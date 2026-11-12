"""
Exact, set-based optimizer for antibody/fluorophore panels.

The objective (spectral matrix complexity) depends only on the selected
fluorophore set, not on which antibody receives each fluorophore. This
implementation therefore enumerates each feasible fluorophore set once,
rather than enumerating every antibody assignment.

Exactness:
- Finds the maximum assignable panel size using bipartite matching.
- Searches every distinct fluorophore set of that size that can be matched
  to distinct antibodies.
- Evaluates every such set and returns the globally best top_n sets.
- No random sampling, heuristic cutoff, or result truncation is used.

The search is still combinatorial in the worst case. This removes duplicate
assignment work and repeated pandas/matrix setup, but cannot make an
exponential optimization problem polynomial-time.
"""

from __future__ import annotations

import csv
import os

# CPU budget: edit this value to cap the number of worker processes.
MAX_WORKERS = 4 #os.cpu_count() or 1

# Avoid BLAS/OpenMP thread oversubscription inside multiple worker processes.
for _thread_var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_thread_var] = "1"

from concurrent.futures import ProcessPoolExecutor
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


def load_antibodies(csv_path: str) -> Dict[str, List[str]]:
    """Read Antibody/Fluorophore pairs from CSV, preserving input order."""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    antibodies: Dict[str, List[str]] = {}
    with open(csv_path, "r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"Antibody", "Fluorophore"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Input CSV must contain Antibody and Fluorophore columns")

        for row in reader:
            antibody = (row["Antibody"] or "").strip()
            fluor = (row["Fluorophore"] or "").strip()
            if not antibody or not fluor:
                continue
            options = antibodies.setdefault(antibody, [])
            if fluor not in options:
                options.append(fluor)
    return antibodies


def normalize_fluorophore(name: str) -> str:
    """Normalize fluorophore names for similarity-matrix matching."""
    return name.strip().replace("/", "-").replace(" ", "").lower()


def load_similarity_matrix(sim_path: str) -> Tuple[pd.DataFrame, List[str]]:
    """Load the similarity matrix; its first column contains row names."""
    df = pd.read_csv(sim_path)
    first_col = df.columns[0]
    df = df.set_index(first_col)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    df = df.apply(pd.to_numeric, errors="coerce")
    return df, df.index.tolist()


def load_detector_map(det_path: str) -> Dict[str, List[str]]:
    """Load detector-to-fluorophore mapping."""
    df = pd.read_csv(det_path)
    result: Dict[str, List[str]] = {}
    for _, row in df.iterrows():
        detector = str(row["Detector"]).strip()
        fluor = str(row["Fluorophore"]).strip()
        result.setdefault(detector, []).append(fluor)
    return result


def fluor_to_similarity_name(
    name: str,
    sim_names_normalized: Dict[str, str],
    det_map: Optional[Dict[str, List[str]]] = None,
) -> Tuple[Optional[str], bool]:
    """Return (matrix name, is_detector_substitution)."""
    normalized = normalize_fluorophore(name)
    if normalized in sim_names_normalized:
        return sim_names_normalized[normalized], False

    if det_map:
        for fluor_list in det_map.values():
            if any(normalize_fluorophore(f) == normalized for f in fluor_list):
                for candidate in fluor_list:
                    candidate_norm = normalize_fluorophore(candidate)
                    if candidate_norm in sim_names_normalized:
                        return sim_names_normalized[candidate_norm], True

    return None, False


def _maximum_matching(
    fluorophores: Sequence[str],
    antibodies_by_fluor: Dict[str, Sequence[str]],
) -> Dict[str, str]:
    """
    Match as many supplied fluorophores as possible to distinct antibodies.

    Return mapping fluorophore -> antibody. A complete result for the supplied
    list means that every fluorophore in the list can be assigned.
    """
    antibody_to_fluor: Dict[str, str] = {}

    # Most constrained fluorophores first makes augmenting paths cheaper.
    ordered = sorted(fluorophores, key=lambda f: len(antibodies_by_fluor[f]))

    def augment(fluor: str, seen_antibodies: set) -> bool:
        for antibody in antibodies_by_fluor[fluor]:
            if antibody in seen_antibodies:
                continue
            seen_antibodies.add(antibody)

            previous = antibody_to_fluor.get(antibody)
            if previous is None or augment(previous, seen_antibodies):
                antibody_to_fluor[antibody] = fluor
                return True
        return False

    for fluor in ordered:
        augment(fluor, set())

    return {fluor: antibody for antibody, fluor in antibody_to_fluor.items()}


def _prepare_similarity(
    sim_df: pd.DataFrame,
    sim_names: Sequence[str],
    fluorophores: Iterable[str],
    det_map: Optional[Dict[str, List[str]]],
):
    """
    Map each candidate fluorophore once and prepare a NumPy similarity array.
    Unknown fluorophores retain the original behaviour: off-diagonal similarity
    to other fluorophores is zero, and the diagonal is one.
    """
    sim_names_normalized = {
        normalize_fluorophore(name): name for name in sim_names
    }
    matrix_names = list(sim_df.index)
    matrix_index = {name: i for i, name in enumerate(matrix_names)}

    # Align rows and columns by name, then symmetrize once.
    aligned = sim_df.reindex(index=matrix_names, columns=matrix_names)
    base = aligned.to_numpy(dtype=np.float64)
    base = (base + base.T) / 2.0

    mapped_index: Dict[str, int] = {}
    substitutions: Dict[str, str] = {}
    unmapped: set = set()

    for fluor in set(fluorophores):
        mapped, substituted = fluor_to_similarity_name(
            fluor, sim_names_normalized, det_map
        )
        if mapped is not None and mapped in matrix_index:
            mapped_index[fluor] = matrix_index[mapped]
            if substituted:
                substitutions[fluor] = mapped
        else:
            mapped_index[fluor] = -1
            unmapped.add(fluor)

    return base, mapped_index, substitutions, unmapped


def _complexity_for_set(
    fluor_set: Sequence[str],
    base_similarity: np.ndarray,
    mapped_index: Dict[str, int],
) -> float:
    """Calculate the same eigenvalue-ratio complexity used in the original."""
    n = len(fluor_set)
    if n < 2:
        return 1.0

    indices = np.fromiter(
        (mapped_index[f] for f in fluor_set), dtype=np.int64, count=n
    )
    known = indices >= 0

    S = np.zeros((n, n), dtype=np.float64)
    if known.any():
        positions = np.flatnonzero(known)
        known_indices = indices[known]
        S[np.ix_(positions, positions)] = base_similarity[
            np.ix_(known_indices, known_indices)
        ]

    np.fill_diagonal(S, 1.0)

    if not np.isfinite(S).all():
        return float("inf")

    eigenvalues = np.linalg.eigvalsh(S)
    smallest = eigenvalues[0]
    if not np.isfinite(eigenvalues).all() or smallest <= 0:
        return float("inf")

    return float(np.sqrt(eigenvalues[-1] / smallest))


_WORKER_DATA = None


def _init_search_worker(data) -> None:
    """Initialize read-only inputs once in each worker process."""
    global _WORKER_DATA
    _WORKER_DATA = data
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(limits=1)
    except ImportError:
        pass


def _search_subtree(task):
    """Search one disjoint prefix subtree; return local top_n and count."""
    start, prefix = task
    d = _WORKER_DATA
    all_fluors, max_count = d["all_fluors"], d["max_count"]
    antibodies_by_fluor, original_order = d["antibodies_by_fluor"], d["original_order"]
    base, mapped_index = d["base"], d["mapped_index"]
    substitutions, unmapped, top_n = d["substitutions"], d["unmapped"], d["top_n"]
    best, evaluated = [], 0

    def search(next_start, chosen):
        nonlocal evaluated
        needed = max_count - len(chosen)
        if needed == 0:
            matching = _maximum_matching(chosen, antibodies_by_fluor)
            if len(matching) != max_count:
                return
            evaluated += 1
            fluor_set = tuple(chosen)
            complexity = _complexity_for_set(fluor_set, base, mapped_index)
            by_antibody = {matching[f]: f for f in matching}
            assignments = {a: by_antibody[a] for a in original_order if a in by_antibody}
            used = set(fluor_set)
            panel = {
                "assignments": assignments,
                "unassigned": [a for a in original_order if a not in assignments],
                "unused_fluorophores": [f for f in all_fluors if f not in used],
                "fluorophores": list(fluor_set),
                "matrix_complexity": complexity,
                "unmapped_fluors": set(f for f in fluor_set if f in unmapped),
                "substitutions": {f: substitutions[f] for f in fluor_set if f in substitutions},
            }
            best.append((complexity, fluor_set, panel))
            best.sort(key=lambda x: (x[0], x[1]))
            if len(best) > top_n:
                best.pop()
            return

        if len(all_fluors) - next_start < needed:
            return
        for i in range(next_start, len(all_fluors) - needed + 1):
            chosen.append(all_fluors[i])
            if len(_maximum_matching(chosen, antibodies_by_fluor)) == len(chosen):
                search(i + 1, chosen)
            chosen.pop()

    search(start, list(prefix))
    return best, evaluated


def _make_search_tasks(all_fluors, antibodies_by_fluor, max_count, target):
    """Create disjoint feasible-prefix subtrees, without listing full panels."""
    frontier = [(0, ())]
    for _ in range(min(max_count, 4)):
        expanded = []
        for start, prefix in frontier:
            needed = max_count - len(prefix)
            if needed <= 0 or len(all_fluors) - start < needed:
                continue
            for i in range(start, len(all_fluors) - needed + 1):
                candidate = prefix + (all_fluors[i],)
                if len(_maximum_matching(candidate, antibodies_by_fluor)) == len(candidate):
                    expanded.append((i + 1, candidate))
        if not expanded:
            return frontier
        frontier = expanded
        if len(frontier) >= target:
            break
    return frontier


def find_best_panels(
    antibodies: Optional[Dict[str, List[str]]] = None,
    csv_path: Optional[str] = None,
    sim_df: Optional[pd.DataFrame] = None,
    sim_names: Optional[Sequence[str]] = None,
    det_map: Optional[Dict[str, List[str]]] = None,
    sim_path: Optional[str] = None,
    det_path: Optional[str] = None,
    top_n: int = 10,
    max_workers: Optional[int] = None,
) -> dict:
    """Find the exact global top_n panels using multiple CPU worker processes."""
    if antibodies is None:
        if not csv_path:
            raise ValueError("Provide antibodies or csv_path")
        antibodies = load_antibodies(csv_path)
    if top_n < 1:
        raise ValueError("top_n must be at least 1")

    workers = max(1, int(MAX_WORKERS if max_workers is None else max_workers))
    options = {a: list(dict.fromkeys(fs)) for a, fs in antibodies.items()}
    original_order = list(options)
    antibodies_by_fluor = {}
    for antibody, fluors in options.items():
        for fluor in fluors:
            antibodies_by_fluor.setdefault(fluor, []).append(antibody)

    all_fluors = sorted(antibodies_by_fluor)
    max_count = len(_maximum_matching(all_fluors, antibodies_by_fluor))
    if max_count == 0:
        return {
            "max_assigned": 0, "total_antibodies": len(options),
            "num_combinations": 1, "num_unique_fluor_sets": 1,
            "combinations": [{"assignments": {}, "unassigned": original_order[:],
                              "unused_fluorophores": all_fluors[:]}],
            "truncated": False,
        }

    if sim_df is None:
        if not sim_path:
            raise ValueError("Provide sim_df or sim_path")
        sim_df, loaded_names = load_similarity_matrix(sim_path)
        if sim_names is None:
            sim_names = loaded_names
    elif sim_names is None:
        sim_names = list(sim_df.index)
    if det_map is None and det_path:
        det_map = load_detector_map(det_path)

    base, mapped_index, substitutions, unmapped = _prepare_similarity(
        sim_df, sim_names, all_fluors, det_map
    )
    data = {
        "all_fluors": all_fluors, "max_count": max_count,
        "antibodies_by_fluor": antibodies_by_fluor,
        "original_order": original_order, "base": base,
        "mapped_index": mapped_index, "substitutions": substitutions,
        "unmapped": unmapped, "top_n": top_n,
    }

    # Multiple moderately sized tasks per worker helps balance uneven branches.
    tasks = _make_search_tasks(all_fluors, antibodies_by_fluor, max_count, workers * 4)
    if workers == 1 or len(tasks) <= 1:
        _init_search_worker(data)
        results = [_search_subtree(t) for t in tasks]
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=_init_search_worker, initargs=(data,)
        ) as pool:
            results = list(pool.map(_search_subtree, tasks, chunksize=1))

    evaluated = sum(count for _, count in results)
    best = [item for local, _ in results for item in local]
    best.sort(key=lambda x: (x[0], x[1]))
    best = best[:top_n]
    panels = [item[2] for item in best]
    combinations = [{
        "assignments": panel["assignments"],
        "unassigned": panel["unassigned"],
        "unused_fluorophores": panel["unused_fluorophores"],
    } for panel in panels]

    return {
        "max_assigned": max_count,
        "total_antibodies": len(options),
        "num_combinations": evaluated,
        "num_unique_fluor_sets": evaluated,
        "combinations": combinations,
        "ranked_panels": panels,
        "truncated": False,
    }


if __name__ == "__main__":
    import time
    base_dir = os.path.dirname(__file__)
    input_path = os.path.join(base_dir, "Data", "input.csv")
    sim_path = os.path.join(base_dir, "constants", "Similarity_Index.csv")
    det_path = os.path.join(base_dir, "constants", "Detector_Fluorophores.csv")

    start = time.time()
    
    antibodies = load_antibodies(input_path)
    sim_df, sim_names = load_similarity_matrix(sim_path)
    det_map = load_detector_map(det_path)

    print(f"Using up to {MAX_WORKERS} worker processes")

    result = find_best_panels(
        antibodies=antibodies,
        sim_df=sim_df,
        sim_names=sim_names,
        det_map=det_map,
        top_n=10,
    )

    elapsed = time.time() - start
    print(
        f"Maximum panel size: "
        f"{result['max_assigned']}/{result['total_antibodies']}"
    )
    print(
        f"Unique feasible fluorophore sets evaluated: "
        f"{result['num_unique_fluor_sets']:,}"
    )
    print(f"Time: {elapsed:.2f}s")
    print("\n--- Globally Best Panels (Lowest Complexity) ---")

    for rank, panel in enumerate(result["ranked_panels"], 1):
        print(
            f"\nRank {rank} - Matrix CI: "
            f"{panel['matrix_complexity']:.4f}"
        )
        for antibody, fluor in panel["assignments"].items():
            tags = []
            if fluor in panel["substitutions"]:
                tags.append(f"[substituted: {panel['substitutions'][fluor]}]")
            if fluor in panel["unmapped_fluors"]:
                tags.append("[unmapped]")
            suffix = " " + " ".join(tags) if tags else ""
            print(f"  {antibody:12} -> {fluor}{suffix}")
        for antibody in panel["unassigned"]:
            print(f"  {antibody:12} -> [unassigned]")
