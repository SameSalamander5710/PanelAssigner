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
    Find globally best spectral panels among all maximum-size feasible panels.

    Unlike find_all_maximum_assignments(), this searches unique fluorophore
    sets, not all antibody/fluorophore permutations. `num_combinations` in the
    returned dictionary means the number of distinct feasible fluorophore
    sets evaluated (not the number of antibody assignments).

    Every feasible set of maximum panel size is evaluated, so top_n is exact.
    """
    if antibodies is None:
        if not csv_path:
            raise ValueError("Provide antibodies or csv_path")
        antibodies = load_antibodies(csv_path)

    if top_n < 1:
        raise ValueError("top_n must be at least 1")

    options = {
        antibody: list(dict.fromkeys(fluors))
        for antibody, fluors in antibodies.items()
    }
    original_order = list(options)

    # Reverse graph: each fluorophore -> antibodies that can use it.
    antibodies_by_fluor: Dict[str, List[str]] = {}
    for antibody, fluors in options.items():
        for fluor in fluors:
            antibodies_by_fluor.setdefault(fluor, []).append(antibody)

    all_fluors = sorted(antibodies_by_fluor)

    # Maximum cardinality matching establishes the required panel size.
    maximum_matching = _maximum_matching(all_fluors, antibodies_by_fluor)
    max_count = len(maximum_matching)

    if max_count == 0:
        return {
            "max_assigned": 0,
            "total_antibodies": len(options),
            "num_combinations": 1,
            "num_unique_fluor_sets": 1,
            "combinations": [{
                "assignments": {},
                "unassigned": original_order[:],
                "unused_fluorophores": all_fluors[:],
            }],
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

    # Keep only the best top_n records; never retain all candidate panels.
    best: List[Tuple[float, Tuple[str, ...], dict]] = []
    evaluated = 0

    # DFS over sorted fluorophore indices generates each set exactly once.
    # A prefix is pruned if its fluorophores cannot be matched to distinct
    # antibodies. Such a prefix can never become a feasible full panel.
    def search(start: int, chosen: List[str]) -> None:
        nonlocal evaluated

        needed = max_count - len(chosen)
        if needed == 0:
            matching = _maximum_matching(chosen, antibodies_by_fluor)
            if len(matching) != max_count:
                return

            evaluated += 1
            fluor_set = tuple(chosen)
            complexity = _complexity_for_set(
                fluor_set, base, mapped_index
            )

            assignments = {
                matching[fluor]: fluor for fluor in matching
            }
            assignments = {
                antibody: assignments[antibody]
                for antibody in original_order
                if antibody in assignments
            }
            unassigned = [
                antibody for antibody in original_order
                if antibody not in assignments
            ]
            used = set(fluor_set)
            solution = {
                "assignments": assignments,
                "unassigned": unassigned,
                "unused_fluorophores": [
                    fluor for fluor in all_fluors if fluor not in used
                ],
                "fluorophores": list(fluor_set),
                "matrix_complexity": complexity,
                "unmapped_fluors": set(f for f in fluor_set if f in unmapped),
                "substitutions": {
                    f: substitutions[f] for f in fluor_set if f in substitutions
                },
            }

            best.append((complexity, fluor_set, solution))
            best.sort(key=lambda item: (item[0], item[1]))
            if len(best) > top_n:
                best.pop()
            return

        # Not enough remaining fluorophores to complete the set.
        if len(all_fluors) - start < needed:
            return

        last_start = len(all_fluors) - needed
        for i in range(start, last_start + 1):
            fluor = all_fluors[i]
            chosen.append(fluor)

            # Feasibility of the prefix is necessary for feasibility of any
            # extension. This pruning does not discard any valid optimum.
            if len(_maximum_matching(chosen, antibodies_by_fluor)) == len(chosen):
                search(i + 1, chosen)

            chosen.pop()

    search(0, [])

    combinations = []
    for _, _, panel in best:
        combinations.append({
            "assignments": panel["assignments"],
            "unassigned": panel["unassigned"],
            "unused_fluorophores": panel["unused_fluorophores"],
        })

    return {
        "max_assigned": max_count,
        "total_antibodies": len(options),
        "num_combinations": evaluated,
        "num_unique_fluor_sets": evaluated,
        "combinations": combinations,
        "ranked_panels": [item[2] for item in best],
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
