"""
Exact, set-based optimizer for antibody/fluorophore panels (branch and bound).

The objective (spectral matrix complexity) depends only on the selected
fluorophore set, not on which antibody receives each fluorophore, so each
feasible fluorophore set is scored once.

Exactness and pruning:
- The maximum panel size comes from bipartite matching.
- Every fluorophore set of that size that can be matched to distinct
  antibodies is a candidate; the globally best top_n sets are returned.
- The complexity sqrt(lambda_max / lambda_min) of a set is computed from a
  principal submatrix of one global similarity matrix. By Cauchy interlacing,
  adding a fluorophore can never lower it. A partial set that is already
  worse than the current n-th best therefore cannot lead to a top_n panel,
  and its whole branch is skipped. This is exact: no heuristic cut-off, no
  sampling; ties are resolved exactly as in an exhaustive search.
- A small tolerance (1e-9 relative) keeps floating-point noise from ever
  discarding a branch that could still tie for a top_n place.

`search_subtree` and `make_search_tasks` are public so the same search can be
split across several workers (see the PyScript page's worker.py).
"""

from __future__ import annotations

import csv
import os
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


INF = float("inf")
_RTOL = 1e-9
_ATOL = 1e-12


def _worse_than(complexity: float, limit: float) -> bool:
    """True if `complexity` is provably worse than `limit` (with tolerance)."""
    return limit != INF and complexity > limit * (1.0 + _RTOL) + _ATOL


def prepare_search(
    antibodies: Dict[str, List[str]],
    sim_df: pd.DataFrame,
    sim_names: Optional[Sequence[str]],
    det_map: Optional[Dict[str, List[str]]],
    top_n: int,
) -> Optional[dict]:
    """Build the read-only search data. Returns None if nothing is assignable."""
    if top_n < 1:
        raise ValueError("top_n must be at least 1")
    options = {a: list(dict.fromkeys(fs)) for a, fs in antibodies.items()}
    original_order = list(options)

    antibodies_by_fluor: Dict[str, List[str]] = {}
    for antibody, fluors in options.items():
        for fluor in fluors:
            antibodies_by_fluor.setdefault(fluor, []).append(antibody)
    all_fluors = sorted(antibodies_by_fluor)

    max_count = len(_maximum_matching(all_fluors, antibodies_by_fluor))
    if max_count == 0:
        return None

    if sim_names is None:
        sim_names = list(sim_df.index)
    base, mapped_index, substitutions, unmapped = _prepare_similarity(
        sim_df, sim_names, all_fluors, det_map
    )
    return {
        "all_fluors": all_fluors, "max_count": max_count,
        "antibodies_by_fluor": antibodies_by_fluor,
        "original_order": original_order, "base": base,
        "mapped_index": mapped_index, "substitutions": substitutions,
        "unmapped": unmapped, "top_n": top_n,
    }


def seed_bound(data: dict) -> float:
    """
    Cheap upper bound on the top_n-th best complexity.

    Builds one panel greedily (feasible sets form a matroid, so a greedy panel
    always reaches full size), then scores all single-fluorophore swaps. The
    top_n-th best of those real, feasible panels bounds the true top_n-th
    best. Returns inf if fewer than top_n panels were found.
    """
    afb, fl = data["antibodies_by_fluor"], data["all_fluors"]
    k, top_n = data["max_count"], data["top_n"]
    base, mi = data["base"], data["mapped_index"]

    def feasible(items):
        return len(_maximum_matching(items, afb)) == len(items)

    chosen: List[str] = []
    while len(chosen) < k:
        pick, pick_c = None, INF
        for f in fl:
            if f in chosen or not feasible(chosen + [f]):
                continue
            c = _complexity_for_set(chosen + [f], base, mi)
            if pick is None or c < pick_c:
                pick, pick_c = f, c
        chosen.append(pick)

    found = {frozenset(chosen): _complexity_for_set(chosen, base, mi)}
    for pos in range(k):
        rest = chosen[:pos] + chosen[pos + 1:]
        for f in fl:
            if f in chosen:
                continue
            cand = rest + [f]
            key = frozenset(cand)
            if key not in found and feasible(cand):
                found[key] = _complexity_for_set(cand, base, mi)
    values = sorted(found.values())
    return values[top_n - 1] if len(values) >= top_n else INF


def next_bound(bound: float, known: Iterable[float], top_n: int) -> float:
    """Tighten `bound` using complexities of panels already found."""
    values = sorted(known)
    if len(values) >= top_n and values[top_n - 1] < bound:
        return values[top_n - 1]
    return bound


def make_search_tasks(data: dict, target: int = 32) -> list:
    """Disjoint feasible-prefix subtrees, most promising (lowest prefix score) first."""
    all_fluors, afb, k = data["all_fluors"], data["antibodies_by_fluor"], data["max_count"]
    frontier = [(0, ())]
    for _ in range(min(k, 4)):
        expanded = []
        for start, prefix in frontier:
            needed = k - len(prefix)
            if needed <= 0 or len(all_fluors) - start < needed:
                continue
            for i in range(start, len(all_fluors) - needed + 1):
                candidate = prefix + (all_fluors[i],)
                if len(_maximum_matching(candidate, afb)) == len(candidate):
                    expanded.append((i + 1, candidate))
        if not expanded:
            break
        frontier = expanded
        if len(frontier) >= target:
            break
    base, mi = data["base"], data["mapped_index"]
    return sorted(
        frontier,
        key=lambda t: (_complexity_for_set(t[1], base, mi) if t[1] else 0.0, t[1]),
    )


def search_subtree(data: dict, task: tuple, bound: float = INF):
    """
    Search one prefix subtree. `bound` is a known upper limit on the top_n-th
    best complexity (inf = none). Returns (best, evaluated, pruned) where best
    is a list of (complexity, fluorophore tuple, panel dict), at most top_n.
    """
    start, prefix = task
    all_fluors, k = data["all_fluors"], data["max_count"]
    afb, original_order = data["antibodies_by_fluor"], data["original_order"]
    base, mi = data["base"], data["mapped_index"]
    substitutions, unmapped, top_n = data["substitutions"], data["unmapped"], data["top_n"]
    best: list = []
    evaluated = 0
    pruned = 0

    def limit() -> float:
        return min(bound, best[-1][0]) if len(best) >= top_n else bound

    def leaf(chosen, complexity, matching) -> None:
        nonlocal evaluated
        evaluated += 1
        fluor_set = tuple(chosen)
        if len(best) >= top_n and (complexity, fluor_set) >= (best[-1][0], best[-1][1]):
            return
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
        best.sort(key=lambda item: (item[0], item[1]))
        if len(best) > top_n:
            best.pop()

    def search(next_start, chosen) -> None:
        nonlocal pruned
        needed = k - len(chosen)
        if len(all_fluors) - next_start < needed:
            return
        for i in range(next_start, len(all_fluors) - needed + 1):
            chosen.append(all_fluors[i])
            matching = _maximum_matching(chosen, afb)
            if len(matching) == len(chosen):
                complexity = _complexity_for_set(chosen, base, mi)
                if _worse_than(complexity, limit()):
                    pruned += 1  # complexity never decreases when adding more
                elif needed == 1:
                    leaf(chosen, complexity, matching)
                else:
                    search(i + 1, chosen)
            chosen.pop()

    chosen = list(prefix)
    if not chosen:
        search(start, chosen)
    else:
        matching = _maximum_matching(chosen, afb)
        if len(matching) == len(chosen):
            complexity = _complexity_for_set(chosen, base, mi)
            if _worse_than(complexity, limit()):
                pruned += 1
            elif len(chosen) == k:
                leaf(chosen, complexity, matching)
            else:
                search(start, chosen)
    return best, evaluated, pruned


def find_best_panels(
    antibodies: Optional[Dict[str, List[str]]] = None,
    csv_path: Optional[str] = None,
    sim_df: Optional[pd.DataFrame] = None,
    sim_names: Optional[Sequence[str]] = None,
    det_map: Optional[Dict[str, List[str]]] = None,
    sim_path: Optional[str] = None,
    det_path: Optional[str] = None,
    top_n: int = 3,
) -> dict:
    """
    Find the exact global top_n panels among all maximum-size feasible panels.

    `num_unique_fluor_sets` counts the full sets actually scored; branches
    proven worse than the current top_n are skipped (`num_pruned_branches`).
    """
    if antibodies is None:
        if not csv_path:
            raise ValueError("Provide antibodies or csv_path")
        antibodies = load_antibodies(csv_path)
    if top_n < 1:
        raise ValueError("top_n must be at least 1")

    if sim_df is None:
        if not sim_path:
            raise ValueError("Provide sim_df or sim_path")
        sim_df, loaded_names = load_similarity_matrix(sim_path)
        if sim_names is None:
            sim_names = loaded_names
    if det_map is None and det_path:
        det_map = load_detector_map(det_path)

    data = prepare_search(antibodies, sim_df, sim_names, det_map, top_n)
    if data is None:
        original_order = list(antibodies)
        return {
            "max_assigned": 0, "total_antibodies": len(antibodies),
            "num_combinations": 1, "num_unique_fluor_sets": 1,
            "num_pruned_branches": 0,
            "combinations": [{
                "assignments": {}, "unassigned": original_order[:],
                "unused_fluorophores": sorted({f for fs in antibodies.values() for f in fs}),
            }],
            "truncated": False,
        }

    bound = seed_bound(data)
    results, known = [], []
    evaluated = pruned = 0
    for task in make_search_tasks(data):
        best, ev, pr = search_subtree(data, task, bound)
        results.extend(best)
        evaluated += ev
        pruned += pr
        known.extend(item[0] for item in best)
        bound = next_bound(bound, known, top_n)

    results.sort(key=lambda item: (item[0], item[1]))
    panels = [item[2] for item in results[:top_n]]
    return {
        "max_assigned": data["max_count"],
        "total_antibodies": len(data["original_order"]),
        "num_combinations": evaluated,
        "num_unique_fluor_sets": evaluated,
        "num_pruned_branches": pruned,
        "combinations": [{
            "assignments": p["assignments"], "unassigned": p["unassigned"],
            "unused_fluorophores": p["unused_fluorophores"],
        } for p in panels],
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

    result = find_best_panels(
        antibodies=antibodies,
        sim_df=sim_df,
        sim_names=sim_names,
        det_map=det_map,
        top_n=3,
    )

    elapsed = time.time() - start
    print(
        f"Maximum panel size: "
        f"{result['max_assigned']}/{result['total_antibodies']}"
    )
    print(
        f"Fluorophore sets scored: {result['num_unique_fluor_sets']:,} "
        f"(branches pruned: {result['num_pruned_branches']:,})"
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
