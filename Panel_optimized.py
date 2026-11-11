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
import re
from functools import lru_cache
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


_SEPARATORS = r"[\s\-/]"
_NUMBER = r"\d+(?:\.\d+)?[A-Za-z]?"  # 421, 5.5, 575V, 647N
_TOKEN_START = r"(?<![A-Za-z0-9])"
_TOKEN_END = r"(?![A-Za-z0-9])"
_TRADEMARK_MARKS = re.compile("[\u00ae\u2122\u00a9]")  # (R), (TM), (C) in vendor names


def load_fluorophore_aliases(csv_path: str) -> List[Tuple[str, str]]:
    """Read Alias,Canonical rows; the header is optional.

    Canonical is the common form every spelling is rewritten to. The input,
    the similarity matrix and the detector map all go through the same rules
    (see `prepare_search`), so the spelling used in the reference files does
    not matter. A trailing `{n}` stands for the number, e.g.
    `AF {n}` -> `Alexa Fluor {n}`.
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Alias CSV not found: {csv_path}")

    rules: List[Tuple[str, str]] = []
    with open(csv_path, "r", newline="", encoding="utf-8-sig") as handle:
        for line_no, row in enumerate(csv.reader(handle), 1):
            if len(row) < 2:
                continue
            alias, canonical = row[0].strip(), row[1].strip()
            if not alias or not canonical:
                continue
            if alias.lower() == "alias" and canonical.lower() == "canonical":
                continue
            problem = _alias_rule_problem(alias, canonical)
            if problem:
                raise ValueError(f"{csv_path}, line {line_no}: {problem}")
            rules.append((alias, canonical))
    return rules


def _alias_rule_problem(alias: str, canonical: str) -> Optional[str]:
    """Describe what is wrong with an alias rule, or None if it is fine."""
    if alias.count("{n}") > 1 or canonical.count("{n}") > 1:
        return "{n} may appear only once"
    if ("{n}" in alias) != ("{n}" in canonical):
        return "{n} must appear in both Alias and Canonical, or in neither"
    if "{n}" in alias and not alias.endswith("{n}"):
        return "{n} must be at the end of Alias"
    if not re.sub(_SEPARATORS + "+", "", alias.replace("{n}", "")):
        return "Alias needs some text besides {n}"
    return None


@lru_cache(maxsize=8)
def _compile_aliases(rules: Tuple[Tuple[str, str], ...]):
    """Build one regex for all rules. Returns (regex, {group index: canonical})."""
    unique: Dict[str, Tuple[str, str]] = {}
    for alias, canonical in rules:  # first rule wins for a repeated alias
        key = re.sub(_SEPARATORS + "+", "", alias).lower()
        unique.setdefault(key, (alias, canonical))
    # Longest alias first, so the most specific rule wins at a given position.
    ordered = sorted(unique.items(), key=lambda item: -len(item[0]))

    alternatives, templates = [], {}
    for i, (_, (alias, canonical)) in enumerate(ordered):
        has_number = alias.endswith("{n}")
        words = re.split(_SEPARATORS + "+", (alias[:-3] if has_number else alias).strip())
        pattern = (_SEPARATORS + "*").join(re.escape(w) for w in words if w)
        if has_number:
            pattern += r"[\s\-]*" + f"(?P<n{i}>{_NUMBER})"
        alternatives.append(f"(?P<r{i}>{_TOKEN_START}{pattern}{_TOKEN_END})")
        templates[i] = canonical
    if not alternatives:
        return None, templates
    return re.compile("|".join(alternatives), re.IGNORECASE), templates


def canonicalize_fluorophore(
    name: str, aliases: Optional[Sequence[Tuple[str, str]]] = None
) -> str:
    """Rewrite known alternative spellings to the canonical fluorophore name.

    Vendor trademark marks are always stripped. Alias rules match whole
    tokens, ignoring case and spaces/hyphens/slashes, anywhere in the name,
    so tandems work too (`PE-AF700`, `PE/Cyanine7`). Unknown names pass
    through unchanged.
    """
    cleaned = _TRADEMARK_MARKS.sub("", name).strip()
    if not aliases:
        return cleaned
    regex, templates = _compile_aliases(tuple(map(tuple, aliases)))
    if regex is None:
        return cleaned

    def replace(match: "re.Match[str]") -> str:
        i = int(match.lastgroup[1:])
        number = match.groupdict().get(f"n{i}")
        template = templates[i]
        return template.replace("{n}", number.upper()) if number else template

    return regex.sub(replace, cleaned)


def _canonicalize_reference(
    sim_df: pd.DataFrame,
    sim_names: Optional[Sequence[str]],
    det_map: Optional[Dict[str, List[str]]],
    aliases: Optional[Sequence[Tuple[str, str]]],
):
    """Put the reference names through the same alias rules as the input.

    Input and reference then meet at the same canonical name, however either
    one is spelled (e.g. "Brilliant Violet 421" in the similarity matrix).
    """
    mapping = {
        name: canonicalize_fluorophore(name, aliases)
        for name in list(sim_df.index) + list(sim_df.columns)
    }
    first_seen: Dict[str, str] = {}
    for name in sim_df.index:
        other = first_seen.setdefault(mapping[name], name)
        if other != name:
            raise ValueError(
                f"Similarity matrix names {other!r} and {name!r} "
                f"both mean {mapping[name]!r}; keep only one"
            )
    sim_df = sim_df.rename(index=mapping, columns=mapping)
    if sim_names is not None:
        sim_names = [canonicalize_fluorophore(n, aliases) for n in sim_names]
    if det_map:
        det_map = {
            detector: [canonicalize_fluorophore(f, aliases) for f in fluors]
            for detector, fluors in det_map.items()
        }
    return sim_df, sim_names, det_map


def _normalize_antibody(name: str) -> str:
    """Normalize antibody names for backbone clash detection."""
    return " ".join(name.strip().lower().split())


def load_backbone(csv_path: str) -> Dict[str, str]:
    """Read fixed backbone assignments from CSV; the header is optional.

    Each antibody keeps its first assignment. A missing file means no backbone.
    """
    if not os.path.exists(csv_path):
        return {}

    backbone: Dict[str, str] = {}
    seen_antibodies: set = set()
    with open(csv_path, "r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if len(row) < 2:
                continue
            antibody = (row[0] or "").strip()
            fluor = (row[1] or "").strip()
            if not antibody or not fluor:
                continue
            if antibody.lower() == "antibody" and fluor.lower() == "fluorophore":
                continue
            key = _normalize_antibody(antibody)
            if key in seen_antibodies:
                continue
            seen_antibodies.add(key)
            backbone[antibody] = fluor
    return backbone


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


def _set_complexity(fluor_set: Sequence[str], data: dict) -> float:
    """Complexity of the fixed backbone plus a candidate free fluorophore set."""
    backbone = data.get("backbone_fluors") or []
    return _complexity_for_set(
        list(backbone) + list(fluor_set), data["base"], data["mapped_index"]
    )


def prepare_search(
    antibodies: Dict[str, List[str]],
    sim_df: pd.DataFrame,
    sim_names: Optional[Sequence[str]],
    det_map: Optional[Dict[str, List[str]]],
    top_n: int,
    backbone: Optional[Dict[str, str]] = None,
    aliases: Optional[Sequence[Tuple[str, str]]] = None,
) -> Optional[dict]:
    """Build the read-only search data. Returns None if nothing is assignable.

    Fluorophore names (input, backbone, similarity matrix and detector map)
    are first rewritten to their canonical spelling using `aliases`, see
    `canonicalize_fluorophore`.
    """
    if top_n < 1:
        raise ValueError("top_n must be at least 1")

    antibodies = {
        antibody: [canonicalize_fluorophore(f, aliases) for f in fluors]
        for antibody, fluors in antibodies.items()
    }
    sim_df, sim_names, det_map = _canonicalize_reference(
        sim_df, sim_names, det_map, aliases
    )
    fixed: Dict[str, str] = {}
    seen_backbone_antibodies: set = set()
    for antibody, fluor in (backbone or {}).items():
        key = _normalize_antibody(antibody)
        if key in seen_backbone_antibodies:  # "CD3" vs "cd3": first one wins
            continue
        seen_backbone_antibodies.add(key)
        fixed[antibody] = canonicalize_fluorophore(fluor, aliases)
    backbone = fixed
    seen_backbone_fluors: Dict[str, str] = {}
    for antibody, fluor in backbone.items():
        key = normalize_fluorophore(fluor)
        if key in seen_backbone_fluors:
            raise ValueError(
                f"Backbone fluorophore {fluor!r} is assigned to both "
                f"{seen_backbone_fluors[key]!r} and {antibody!r}"
            )
        seen_backbone_fluors[key] = antibody

    backbone_fluors = list(backbone.values())
    backbone_fluor_keys = set(seen_backbone_fluors)
    backbone_antibody_keys = {_normalize_antibody(a) for a in backbone}

    options: Dict[str, List[str]] = {}
    for antibody, fluors in antibodies.items():
        if _normalize_antibody(antibody) in backbone_antibody_keys:
            continue
        options[antibody] = [
            fluor for fluor in dict.fromkeys(fluors)
            if normalize_fluorophore(fluor) not in backbone_fluor_keys
        ]

    free_original_order = list(options)
    original_order = list(backbone) + free_original_order

    antibodies_by_fluor: Dict[str, List[str]] = {}
    for antibody, fluors in options.items():
        for fluor in fluors:
            antibodies_by_fluor.setdefault(fluor, []).append(antibody)
    all_fluors = sorted(antibodies_by_fluor)

    free_max_count = len(_maximum_matching(all_fluors, antibodies_by_fluor))
    max_count = len(backbone) + free_max_count
    if max_count == 0:
        return None

    if sim_names is None:
        sim_names = list(sim_df.index)
    base, mapped_index, substitutions, unmapped = _prepare_similarity(
        sim_df, sim_names, list(all_fluors) + backbone_fluors, det_map
    )
    return {
        "all_fluors": all_fluors, "max_count": max_count,
        "free_max_count": free_max_count,
        "antibodies_by_fluor": antibodies_by_fluor,
        "original_order": original_order,
        "free_original_order": free_original_order,
        "backbone_assignments": backbone,
        "backbone_fluors": backbone_fluors,
        "base": base,
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
    k, top_n = data["free_max_count"], data["top_n"]

    def feasible(items):
        return len(_maximum_matching(items, afb)) == len(items)

    chosen: List[str] = []
    while len(chosen) < k:
        pick, pick_c = None, INF
        for f in fl:
            if f in chosen or not feasible(chosen + [f]):
                continue
            c = _set_complexity(chosen + [f], data)
            if pick is None or c < pick_c:
                pick, pick_c = f, c
        chosen.append(pick)

    found = {frozenset(chosen): _set_complexity(chosen, data)}
    for pos in range(k):
        rest = chosen[:pos] + chosen[pos + 1:]
        for f in fl:
            if f in chosen:
                continue
            cand = rest + [f]
            key = frozenset(cand)
            if key not in found and feasible(cand):
                found[key] = _set_complexity(cand, data)
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
    all_fluors, afb, k = data["all_fluors"], data["antibodies_by_fluor"], data["free_max_count"]
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
    return sorted(
        frontier,
        key=lambda t: (_set_complexity(t[1], data) if t[1] else 0.0, t[1]),
    )


def search_subtree(data: dict, task: tuple, bound: float = INF):
    """
    Search one prefix subtree. `bound` is a known upper limit on the top_n-th
    best complexity (inf = none). Returns (best, evaluated, pruned) where best
    is a list of (complexity, fluorophore tuple, panel dict), at most top_n.
    """
    start, prefix = task
    all_fluors, k = data["all_fluors"], data["free_max_count"]
    afb = data["antibodies_by_fluor"]
    original_order = data["original_order"]
    free_original_order = data["free_original_order"]
    backbone_assignments = data["backbone_assignments"]
    backbone_fluors = data["backbone_fluors"]
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
        assignments = {}
        for antibody in original_order:
            if antibody in backbone_assignments:
                assignments[antibody] = backbone_assignments[antibody]
            elif antibody in by_antibody:
                assignments[antibody] = by_antibody[antibody]
        full_fluors = list(backbone_fluors) + list(fluor_set)
        used = set(fluor_set)
        panel = {
            "assignments": assignments,
            "unassigned": [a for a in free_original_order if a not in by_antibody],
            "unused_fluorophores": [f for f in all_fluors if f not in used],
            "fluorophores": full_fluors,
            "matrix_complexity": complexity,
            "unmapped_fluors": set(f for f in full_fluors if f in unmapped),
            "substitutions": {f: substitutions[f] for f in full_fluors if f in substitutions},
            "backbone": dict(backbone_assignments),
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
                complexity = _set_complexity(chosen, data)
                if _worse_than(complexity, limit()):
                    pruned += 1  # complexity never decreases when adding more
                elif needed == 1:
                    leaf(chosen, complexity, matching)
                else:
                    search(i + 1, chosen)
            chosen.pop()

    if k == 0:
        leaf([], _set_complexity([], data), {})
        return best, evaluated, pruned

    chosen = list(prefix)
    if not chosen:
        search(start, chosen)
    else:
        matching = _maximum_matching(chosen, afb)
        if len(matching) == len(chosen):
            complexity = _set_complexity(chosen, data)
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
    backbone: Optional[Dict[str, str]] = None,
    aliases: Optional[Sequence[Tuple[str, str]]] = None,
    alias_path: Optional[str] = None,
) -> dict:
    """
    Find the exact global top_n panels among all maximum-size feasible panels.

    `backbone` is a fixed antibody -> fluorophore mapping that is always part
    of every panel and takes priority over conflicting input antibodies or
    fluorophores. `aliases` (or `alias_path`) rewrites alternative fluorophore
    spellings such as "Brilliant Violet 421" to the canonical "BV421" before
    matching. `num_unique_fluor_sets` counts the full sets actually scored;
    branches proven worse than the current top_n are skipped
    (`num_pruned_branches`).
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
    if aliases is None and alias_path:
        aliases = load_fluorophore_aliases(alias_path)

    data = prepare_search(
        antibodies, sim_df, sim_names, det_map, top_n, backbone, aliases
    )
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
    backbone_path = os.path.join(base_dir, "Data", "input_backbone.csv")
    sim_path = os.path.join(base_dir, "constants", "Similarity_Index.csv")
    det_path = os.path.join(base_dir, "constants", "Detector_Fluorophores.csv")
    alias_path = os.path.join(base_dir, "constants", "Fluorophore_Aliases.csv")

    start = time.time()
    
    antibodies = load_antibodies(input_path)
    backbone = load_backbone(backbone_path)
    sim_df, sim_names = load_similarity_matrix(sim_path)
    det_map = load_detector_map(det_path)
    aliases = load_fluorophore_aliases(alias_path)

    result = find_best_panels(
        antibodies=antibodies,
        sim_df=sim_df,
        sim_names=sim_names,
        det_map=det_map,
        top_n=3,
        backbone=backbone,
        aliases=aliases,
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
            if antibody in panel.get("backbone", {}):
                tags.append("[backbone]")
            if fluor in panel["substitutions"]:
                tags.append(f"[substituted: {panel['substitutions'][fluor]}]")
            if fluor in panel["unmapped_fluors"]:
                tags.append("[unmapped]")
            suffix = " " + " ".join(tags) if tags else ""
            print(f"  {antibody:12} -> {fluor}{suffix}")
        for antibody in panel["unassigned"]:
            print(f"  {antibody:12} -> [unassigned]")
