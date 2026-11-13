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
import math
import os
import re
from functools import lru_cache
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


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
    sim_df: SimilarityMatrix,
    sim_names: Optional[Sequence[str]],
    det_map: Optional[Dict[str, List[str]]],
    aliases: Optional[Sequence[Tuple[str, str]]],
):
    """Put the reference names through the same alias rules as the input.

    Input and reference then meet at the same canonical name, however either
    one is spelled (e.g. "Brilliant Violet 421" in the similarity matrix).
    """
    sim_df = _as_matrix(sim_df)
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
    sim_df = sim_df.renamed(mapping)
    if sim_names is not None:
        sim_names = [canonicalize_fluorophore(n, aliases) for n in sim_names]
    if det_map:
        det_map = {
            detector: [canonicalize_fluorophore(f, aliases) for f in fluors]
            for detector, fluors in det_map.items()
        }
    return sim_df, sim_names, det_map


def canonical_types(
    types: Optional[Dict[str, str]],
    aliases: Optional[Sequence[Tuple[str, str]]],
) -> Dict[str, str]:
    """{normalized canonical fluorophore: lower-case Type}, see `prepare_reference`."""
    return {
        normalize_fluorophore(canonicalize_fluorophore(f, aliases)): kind.lower()
        for f, kind in (types or {}).items()
    }


def prepare_reference(
    sim_df: SimilarityMatrix,
    sim_names: Optional[Sequence[str]],
    det_map: Optional[Dict[str, List[str]]],
    aliases: Optional[Sequence[Tuple[str, str]]],
    types: Optional[Dict[str, str]] = None,
):
    """Put the reference data through the alias rules once.

    Returns (sim_df, sim_names, det_map, types) to pass to `prepare_search`,
    `suggest_additions` and `eligible_fluorophore_count` with
    `canonical=True`. Those functions otherwise redo this on every call, which
    is slow when it is repeated, e.g. for each result panel.
    """
    sim_df, sim_names, det_map = _canonicalize_reference(sim_df, sim_names, det_map, aliases)
    return sim_df, sim_names, det_map, canonical_types(types, aliases)


def _normalize_antibody(name: str) -> str:
    """Normalize antibody names for backbone clash detection."""
    return " ".join(name.strip().lower().split())


def load_backbone(csv_path: str) -> Dict[str, str]:
    """Read fixed backbone assignments from CSV; the header is optional.

    Each antibody keeps its first assignment. A row with a single value is a
    bare fluorophore and gets an automatic antibody name (Avoid_1, Avoid_2,
    ...). A missing file means no backbone.
    """
    if not os.path.exists(csv_path):
        return {}

    backbone: Dict[str, str] = {}
    seen_antibodies: set = set()
    auto = 0
    with open(csv_path, "r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        for row in reader:
            fields = [field.strip() for field in row if field.strip()]
            if not fields:
                continue
            if len(fields) >= 2:
                antibody, fluor = fields[0], fields[1]
                if antibody.lower() == "antibody" and fluor.lower() == "fluorophore":
                    continue
            else:
                auto += 1
                antibody, fluor = f"Avoid_{auto}", fields[0]
            key = _normalize_antibody(antibody)
            if key in seen_antibodies:
                continue
            seen_antibodies.add(key)
            backbone[antibody] = fluor
    return backbone


class SimilarityMatrix:
    """A similarity table: row names, column names and a float array.

    This is all the search needs from the matrix. It replaces a pandas
    DataFrame so that the web version only has to load numpy; a DataFrame
    passed to the public functions still works (see `_as_matrix`).
    """

    def __init__(self, index: Sequence[str], columns: Sequence[str], values) -> None:
        self.index: List[str] = list(index)
        self.columns: List[str] = list(columns)
        self.values = np.asarray(values, dtype=np.float64).reshape(len(self.index), len(self.columns))
        self._symmetric: Optional[np.ndarray] = None

    def to_numpy(self) -> np.ndarray:
        return self.values

    def renamed(self, mapping: Dict[str, str]) -> "SimilarityMatrix":
        return SimilarityMatrix(
            [mapping.get(n, n) for n in self.index],
            [mapping.get(n, n) for n in self.columns],
            self.values,
        )

    def aligned(self) -> np.ndarray:
        """Values with the columns put in row order (NaN where a row has no column)."""
        position: Dict[str, int] = {}
        for j, name in enumerate(self.columns):
            position.setdefault(name, j)
        out = np.full((len(self.index), len(self.index)), np.nan)
        for i, name in enumerate(self.index):
            j = position.get(name)
            if j is not None:
                out[:, i] = self.values[:, j]
        return out


    def symmetrized(self) -> np.ndarray:
        """Aligned values averaged with their transpose (computed once, read-only)."""
        if self._symmetric is None:
            aligned = self.aligned()
            self._symmetric = (aligned + aligned.T) / 2.0
            self._symmetric.flags.writeable = False
        return self._symmetric


def _as_matrix(sim) -> SimilarityMatrix:
    """Accept a SimilarityMatrix or anything shaped like a DataFrame."""
    if isinstance(sim, SimilarityMatrix):
        return sim
    return SimilarityMatrix(
        [str(n) for n in sim.index], [str(n) for n in sim.columns], np.asarray(sim, dtype=np.float64)
    )


def _read_table(path: str) -> List[List[str]]:
    """All non-blank CSV rows (a UTF-8 BOM is ignored)."""
    with open(path, "r", newline="", encoding="utf-8-sig") as handle:
        return [row for row in csv.reader(handle) if any(cell.strip() for cell in row)]


def _to_float(text: str) -> float:
    try:
        return float(text)
    except ValueError:
        return float("nan")  # unreadable cells count as missing, as before


def load_similarity_matrix(sim_path: str) -> Tuple[SimilarityMatrix, List[str]]:
    """Load the similarity matrix; its first column contains row names."""
    rows = _read_table(sim_path)
    columns = [cell.strip() for cell in rows[0][1:]]
    index: List[str] = []
    values = np.full((len(rows) - 1, len(columns)), np.nan)
    for i, row in enumerate(rows[1:]):
        index.append(row[0].strip())
        for j, cell in enumerate(row[1:1 + len(columns)]):
            values[i, j] = _to_float(cell)
    return SimilarityMatrix(index, columns, values), list(index)


def _detector_rows(det_path: str):
    """Yield (detector, fluorophore, row dict by column name) for each filled row."""
    rows = _read_table(det_path)
    header = [cell.strip() for cell in rows[0]]
    for row in rows[1:]:
        record = {name: (row[k].strip() if k < len(row) else "") for k, name in enumerate(header)}
        if record.get("Detector") and record.get("Fluorophore"):
            yield record["Detector"], record["Fluorophore"], record


def load_detector_map(det_path: str) -> Dict[str, List[str]]:
    """Load detector-to-fluorophore mapping."""
    result: Dict[str, List[str]] = {}
    for detector, fluor, _ in _detector_rows(det_path):
        result.setdefault(detector, []).append(fluor)
    return result


# Detector_Fluorophores.csv "Type" values that may be suggested as additions.
# Other types (Dye, Fluorescent Protein, Microscopy, ...) are kept in the file
# for substitution and naming but are not offered as antibody labels.
SUGGESTED_TYPES: Tuple[str, ...] = ("Antibody",)


def load_fluorophore_types(det_path: str) -> Dict[str, str]:
    """Load {fluorophore: Type} from the detector CSV's optional Type column.

    Returns an empty dict when the file has no Type column, which means no
    fluorophore is filtered out of the suggestions.
    """
    types: Dict[str, str] = {}
    for _, fluor, record in _detector_rows(det_path):
        if record.get("Type"):
            types.setdefault(fluor, record["Type"])
    return types


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
    sim_df: SimilarityMatrix,
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
    sim_df = _as_matrix(sim_df)
    matrix_names = list(sim_df.index)
    matrix_index = {name: i for i, name in enumerate(matrix_names)}

    # Rows and columns are aligned by name and symmetrized once per matrix.
    base = sim_df.symmetrized()

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


def _similarity_submatrix(
    fluor_set: Sequence[str],
    base_similarity: np.ndarray,
    mapped_index: Dict[str, int],
) -> np.ndarray:
    """Similarity matrix of a fluorophore set, as used for its complexity.

    Unmapped fluorophores have zero similarity to all others; the diagonal
    is always one.
    """
    n = len(fluor_set)
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
    return S


def panel_similarity_matrix(data: dict, fluor_set: Sequence[str]) -> np.ndarray:
    """The exact matrix behind a panel's complexity score, rows/columns in
    the order of `fluor_set` (`data` comes from `prepare_search`)."""
    return _similarity_submatrix(fluor_set, data["base"], data["mapped_index"])


def _complexity_for_set(
    fluor_set: Sequence[str],
    base_similarity: np.ndarray,
    mapped_index: Dict[str, int],
) -> float:
    """Calculate the same eigenvalue-ratio complexity used in the original."""
    if len(fluor_set) < 2:
        return 1.0

    S = _similarity_submatrix(fluor_set, base_similarity, mapped_index)

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


_PAIR_SAFE = 1.0 - 1e-9  # closer to 1 than this, leave the pair to the eigenvalue solve


class _SetComplexity:
    """Complexity of a fixed set plus fluorophores added one after another.

    The similarity matrix of the whole set is kept in a work array that grows
    row by row, so adding a fluorophore costs one row lookup and one
    eigenvalue solve, instead of rebuilding the matrix from names. `M` is the
    full similarity table (see `_similarity_submatrix`) of every fluorophore
    that can occur, and fluorophores are given by their row in it.
    """

    def __init__(self, M: np.ndarray, base_rows: Sequence[int], max_extra: int) -> None:
        self.M = M
        self.nb = len(base_rows)
        size = self.nb + max_extra
        self.rows = np.zeros(size, dtype=np.intp)
        self.rows[: self.nb] = base_rows
        self.W = np.zeros((size, size))  # only the lower triangle is used
        for r in range(self.nb):
            self.W[r, :r] = M[base_rows[r], self.rows[:r]]
            self.W[r, r] = 1.0
        self.finite = bool(np.isfinite(M).all())
        self.signed = bool((M < 0).any())

    def push(self, depth: int, row: int) -> float:
        """Make `row` the added fluorophore number `depth` (0 = first).

        Returns its largest similarity to the fluorophores before it. By
        eigenvalue interlacing, a set containing a pair with similarity s has
        complexity of at least sqrt((1 + s) / (1 - s)), so this can rule a
        set out without solving for eigenvalues.
        """
        m = self.nb + depth
        self.rows[m] = row
        self.W[m, m] = 1.0
        if m == 0:
            return 0.0
        values = self.M[row, self.rows[:m]]
        self.W[m, :m] = values
        return float(np.abs(values).max() if self.signed else values.max())

    def complexity(self, depth: int) -> float:
        """Complexity with the first `depth + 1` added fluorophores (-1: none)."""
        m = self.nb + depth + 1
        if m < 2:
            return 1.0
        S = self.W[:m, :m]
        if not self.finite and not np.isfinite(S).all():
            return INF
        eigenvalues = np.linalg.eigvalsh(S)
        smallest = eigenvalues[0]
        if not np.isfinite(eigenvalues).all() or smallest <= 0:
            return INF
        return float(np.sqrt(eigenvalues[-1] / smallest))


def _pair_limit(limit: float) -> float:
    """Pair similarity above which a set is provably worse than `limit`."""
    if limit == INF:
        return 2.0  # nothing is that similar
    bound = limit * (1.0 + _RTOL) + _ATOL
    return (bound * bound - 1.0) / (bound * bound + 1.0)


def prepare_search(
    antibodies: Dict[str, List[str]],
    sim_df: SimilarityMatrix,
    sim_names: Optional[Sequence[str]],
    det_map: Optional[Dict[str, List[str]]],
    top_n: int,
    backbone: Optional[Dict[str, str]] = None,
    aliases: Optional[Sequence[Tuple[str, str]]] = None,
    max_assigned: Optional[int] = None,
    canonical: bool = False,
) -> Optional[dict]:
    """Build the read-only search data. Returns None if nothing is assignable.

    `canonical=True` says that sim_df, sim_names and det_map already come from
    `prepare_reference`, so they are not rewritten again.

    `max_assigned` caps the panel size (backbone included): the search then
    looks for the best panels fitting at most that many antibodies. None, or
    anything above what is possible, means as many as can be assigned.

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
    if not canonical:
        sim_df, sim_names, det_map = _canonicalize_reference(
            sim_df, sim_names, det_map, aliases
        )
    fixed: Dict[str, str] = {}
    seen_backbone_antibodies: set = set()
    seen_backbone_fluors: Dict[str, str] = {}
    for antibody, fluor in (backbone or {}).items():
        key = _normalize_antibody(antibody)
        if key in seen_backbone_antibodies:  # "CD3" vs "cd3": first one wins
            continue
        fluor = canonicalize_fluorophore(fluor, aliases)
        fluor_key = normalize_fluorophore(fluor)
        if fluor_key in seen_backbone_fluors:  # same fluorophore again: first one wins
            continue
        seen_backbone_antibodies.add(key)
        seen_backbone_fluors[fluor_key] = antibody
        fixed[antibody] = fluor
    backbone = fixed

    backbone_fluors = list(backbone.values())
    backbone_fluor_keys = set(seen_backbone_fluors)
    backbone_antibody_keys = {_normalize_antibody(a) for a in backbone}

    # Treat spellings that normalize to the same fluorophore (e.g. "PE/Cy7"
    # and "PE-Cy7") as one option, using a single spelling everywhere so the
    # same fluorophore can never appear twice or form a duplicate panel.
    canonical_spelling: Dict[str, str] = {}
    options: Dict[str, List[str]] = {}
    for antibody, fluors in antibodies.items():
        if _normalize_antibody(antibody) in backbone_antibody_keys:
            continue
        deduped: List[str] = []
        seen_keys: set = set()
        for fluor in fluors:
            key = normalize_fluorophore(fluor)
            if key in backbone_fluor_keys or key in seen_keys:
                continue
            seen_keys.add(key)
            deduped.append(canonical_spelling.setdefault(key, fluor))
        options[antibody] = deduped

    free_original_order = list(options)
    original_order = list(backbone) + free_original_order

    antibodies_by_fluor: Dict[str, List[str]] = {}
    for antibody, fluors in options.items():
        for fluor in fluors:
            antibodies_by_fluor.setdefault(fluor, []).append(antibody)
    all_fluors = sorted(antibodies_by_fluor)

    free_max_count = len(_maximum_matching(all_fluors, antibodies_by_fluor))
    max_possible = len(backbone) + free_max_count
    if max_assigned is not None:
        if max_assigned < 1:
            raise ValueError('"Fit at most" (max_assigned) must be at least 1.')
        if max_assigned < len(backbone):
            raise ValueError(
                f'"Fit at most" is {max_assigned}, but the fixed backbone '
                f"alone has {len(backbone)} antibodies."
            )
        free_max_count = min(free_max_count, max_assigned - len(backbone))
    max_count = len(backbone) + free_max_count
    if max_count == 0:
        return None

    if sim_names is None:
        sim_names = list(sim_df.index)
    base, mapped_index, substitutions, unmapped = _prepare_similarity(
        sim_df, sim_names, list(all_fluors) + backbone_fluors, det_map
    )
    # Search tables. The search tries the most-overlapping fluorophores first:
    # their conflicts show up near the top of the tree, so more of it is cut
    # off. Results do not depend on this order. `all_fluors` stays
    # alphabetical for everything that is shown.
    labels = list(backbone_fluors) + list(all_fluors)
    M = _similarity_submatrix(labels, base, mapped_index)
    n_free = len(all_fluors)
    overlap = (M[len(backbone_fluors):, len(backbone_fluors):].sum(axis=1) - 1.0) / max(n_free - 1, 1)
    overlap = np.nan_to_num(overlap, nan=0.0)
    ranked = sorted(range(n_free), key=lambda i: (-overlap[i], all_fluors[i]))
    search_order = [all_fluors[i] for i in ranked]
    return {
        "search_order": search_order,
        "search_pos": {f: i for i, f in enumerate(search_order)},
        "label_pos": {f: i for i, f in enumerate(labels)},
        "M": M,
        "all_fluors": all_fluors, "max_count": max_count,
        "max_possible": max_possible,
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
    """Disjoint feasible-prefix subtrees, most promising (lowest prefix score) first.

    A task is (start, prefix): the fluorophores of `prefix` are fixed and the
    rest are taken from `search_order[start:]`.
    """
    order, afb, k = data["search_order"], data["antibodies_by_fluor"], data["free_max_count"]
    frontier = [(0, ())]
    for _ in range(min(k, 4)):
        expanded = []
        for start, prefix in frontier:
            needed = k - len(prefix)
            if needed <= 0 or len(order) - start < needed:
                continue
            for i in range(start, len(order) - needed + 1):
                candidate = prefix + (order[i],)
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

    Fluorophores are tried in `search_order`, but every reported panel, its
    complexity and its tie-breaking use the alphabetical order, so the result
    does not depend on the search order.
    """
    start, prefix = task
    order, k = data["search_order"], data["free_max_count"]
    afb = data["antibodies_by_fluor"]
    all_fluors = data["all_fluors"]
    original_order = data["original_order"]
    free_original_order = data["free_original_order"]
    backbone_assignments = data["backbone_assignments"]
    backbone_fluors = data["backbone_fluors"]
    substitutions, unmapped, top_n = data["substitutions"], data["unmapped"], data["top_n"]
    best: list = []
    evaluated = 0
    pruned = 0

    label_pos = data["label_pos"]
    n_order = len(order)
    rows = [label_pos[f] for f in order]  # row of each fluorophore in M
    options = [afb[f] for f in order]  # its antibodies
    builder = _SetComplexity(data["M"], [label_pos[f] for f in backbone_fluors], k)

    def limit() -> float:
        return min(bound, best[-1][0]) if len(best) >= top_n else bound

    def leaf(chosen: List[int]) -> None:
        nonlocal evaluated
        evaluated += 1
        fluor_set = tuple(sorted(order[i] for i in chosen))
        complexity = _set_complexity(fluor_set, data)  # exact value, alphabetical order
        if len(best) >= top_n and (complexity, fluor_set) >= (best[-1][0], best[-1][1]):
            return
        matching = _maximum_matching(fluor_set, afb)
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

    def add(matched: Dict[str, int], i: int) -> Optional[Dict[str, int]]:
        """The antibody matching of the chosen fluorophores plus `i`, or None.

        `matched` maps antibody -> fluorophore and covers every chosen
        fluorophore, so one augmenting-path search from the new one suffices.
        """
        new = dict(matched)
        for antibody in options[i]:
            if antibody not in new:  # a free antibody: nothing to rearrange
                new[antibody] = i
                return new
        seen: set = set()

        def augment(j: int) -> bool:
            for antibody in options[j]:
                if antibody in seen:
                    continue
                seen.add(antibody)
                holder = new.get(antibody)
                if holder is None or augment(holder):
                    new[antibody] = j
                    return True
            return False

        return new if augment(i) else None

    def search(next_start: int, chosen: List[int], matched: Dict[str, int]) -> None:
        nonlocal pruned
        depth = len(chosen)
        needed = k - depth
        if n_order - next_start < needed:
            return
        for i in range(next_start, n_order - needed + 1):
            grown = add(matched, i)
            if grown is None:
                continue  # no antibody assignment for these fluorophores
            top_similarity = builder.push(depth, rows[i])
            current = limit()
            if current != INF and _PAIR_SAFE > top_similarity > _pair_limit(current):
                pruned += 1  # one pair is already too similar
                continue
            if _worse_than(builder.complexity(depth), current):
                pruned += 1  # complexity never decreases when adding more
            elif needed == 1:
                chosen.append(i)
                leaf(chosen)
                chosen.pop()
            else:
                chosen.append(i)
                search(i + 1, chosen, grown)
                chosen.pop()

    if k == 0:
        leaf([])
        return best, evaluated, pruned

    chosen: List[int] = []
    matched: Optional[Dict[str, int]] = {}
    for f in prefix:
        i = data["search_pos"][f]
        matched = add(matched, i)
        if matched is None:
            return best, evaluated, pruned
        builder.push(len(chosen), rows[i])
        chosen.append(i)
    if chosen:
        complexity = builder.complexity(len(chosen) - 1)
        if _worse_than(complexity, limit()):
            pruned += 1
        elif len(chosen) == k:
            leaf(chosen)
        else:
            search(start, chosen, matched)
    else:
        search(start, chosen, matched)
    return best, evaluated, pruned


MAX_SUGGESTED = 5  # adding more than this at once gets slow and is rarely useful


def _suggestion_pool(
    sim_df: SimilarityMatrix,
    sim_names: Optional[Sequence[str]],
    det_map: Optional[Dict[str, List[str]]],
    panel_fluors: Sequence[str],
    aliases: Optional[Sequence[Tuple[str, str]]] = None,
    types: Optional[Dict[str, str]] = None,
    allowed_types: Sequence[str] = SUGGESTED_TYPES,
    canonical: bool = False,
) -> dict:
    """Prepare the eligible additions shared by the suggestion functions.

    Returns the canonical panel, its similarity data and base complexity, and
    the eligible pool ordered by best single addition first, with that
    complexity alongside.

    `types` ({fluorophore: Type}, see `load_fluorophore_types`) limits the
    candidates to those whose Type is in `allowed_types`. A candidate with no
    Type listed stays eligible, so a missing entry never hides it silently.
    The panel itself is never filtered.

    `canonical=True` says that the reference data and `types` come from
    `prepare_reference`.
    """
    if not canonical:
        sim_df, sim_names, det_map = _canonicalize_reference(sim_df, sim_names, det_map, aliases)
        types = canonical_types(types, aliases)
    if sim_names is None:
        sim_names = list(sim_df.index)

    panel = [canonicalize_fluorophore(f, aliases) for f in panel_fluors]
    taken = {normalize_fluorophore(f) for f in panel}
    candidates = [
        name for name in dict.fromkeys(sim_names)
        if normalize_fluorophore(name) not in taken
    ]
    if types:
        allowed = {kind.lower() for kind in allowed_types}
        candidates = [
            c for c in candidates
            if normalize_fluorophore(c) not in types  # no Type listed: keep
            or types[normalize_fluorophore(c)] in allowed
        ]
    base, mapped_index, substitutions, _ = _prepare_similarity(
        sim_df, sim_names, panel + candidates, det_map
    )
    # A panel fluorophore that stands in for a matrix fluorophore occupies it.
    taken |= {normalize_fluorophore(s) for s in substitutions.values()}
    candidates = [c for c in candidates if normalize_fluorophore(c) not in taken]

    labels = panel + candidates
    M = _similarity_submatrix(labels, base, mapped_index)
    label_pos = {name: len(panel) + i for i, name in enumerate(candidates)}
    builder = _SetComplexity(M, range(len(panel)), 1)

    base_complexity = builder.complexity(-1)
    if not np.isfinite(base_complexity):
        raise ValueError(
            "Two fluorophores in this panel have identical spectra, so its "
            "complexity is infinite and additions cannot be ranked."
        )

    # Candidates that are infinite on their own (identical to a panel member)
    # can never be part of a finite set. Best single additions first, so good
    # sets are found early and the bound tightens quickly.
    singles = []
    for name in candidates:
        builder.push(0, label_pos[name])
        c = builder.complexity(0)
        if np.isfinite(c):
            singles.append((c, name))
    singles.sort()
    return {
        "panel": panel,
        "base": base,
        "mapped_index": mapped_index,
        "M": M,
        "label_pos": label_pos,
        "base_complexity": base_complexity,
        "pool": [name for _, name in singles],
        "single_complexity": [c for c, _ in singles],
    }


def eligible_fluorophore_count(
    sim_df: SimilarityMatrix,
    sim_names: Optional[Sequence[str]],
    det_map: Optional[Dict[str, List[str]]],
    panel_fluors: Sequence[str],
    aliases: Optional[Sequence[Tuple[str, str]]] = None,
    types: Optional[Dict[str, str]] = None,
    allowed_types: Sequence[str] = SUGGESTED_TYPES,
    canonical: bool = False,
) -> int:
    """How many reference fluorophores could still be added to a panel.

    Returns 0 when the panel itself has infinite complexity, since then no
    addition can be ranked.
    """
    try:
        return len(_suggestion_pool(
            sim_df, sim_names, det_map, panel_fluors, aliases, types, allowed_types,
            canonical,
        )["pool"])
    except ValueError:
        return 0


def suggest_additions(
    sim_df: SimilarityMatrix,
    sim_names: Optional[Sequence[str]],
    det_map: Optional[Dict[str, List[str]]],
    panel_fluors: Sequence[str],
    n: int,
    top_n: int = 5,
    aliases: Optional[Sequence[Tuple[str, str]]] = None,
    types: Optional[Dict[str, str]] = None,
    allowed_types: Sequence[str] = SUGGESTED_TYPES,
    canonical: bool = False,
) -> dict:
    """Best sets of `n` extra fluorophores to add to an existing panel.

    Candidates are all fluorophores of the reference matrix that the panel
    does not use yet (only those of an allowed Type, when `types` is given,
    so dyes and fluorescent proteins are not suggested as antibody labels).
    Sets are ranked by the complexity of the panel plus the
    new fluorophores. This is exact: the complexity never decreases when a
    fluorophore is added, so partial sets that are already worse than the
    current top_n are skipped. (The best n+1 set does not always contain the
    best n set, so adding one at a time would not be exact.)

    Returns {"base_complexity", "candidates", "combinations", "suggestions":
    [{"fluors", "matrix_complexity", "values"}]}, where `candidates` is the
    number of eligible fluorophores, `combinations` is how many `n`-sized sets
    they form, and `values` is the similarity matrix of the panel followed by
    the new fluorophores.
    """
    if not 1 <= n <= MAX_SUGGESTED:
        raise ValueError(f"Choose between 1 and {MAX_SUGGESTED} fluorophores to add.")
    context = _suggestion_pool(
        sim_df, sim_names, det_map, panel_fluors, aliases, types, allowed_types,
        canonical,
    )
    panel = context["panel"]
    base = context["base"]
    mapped_index = context["mapped_index"]
    base_complexity = context["base_complexity"]
    pool = context["pool"]
    single_complexity = context["single_complexity"]

    if len(pool) < n:
        raise ValueError(
            f"Only {len(pool)} suitable fluorophores are left in the reference "
            f"matrix to add."
        )

    found: List[Tuple[float, Tuple[str, ...]]] = []
    builder = _SetComplexity(context["M"], range(len(panel)), n)
    pool_rows = [context["label_pos"][name] for name in pool]

    def limit() -> float:
        return found[-1][0] if len(found) >= top_n else INF

    def extend(start: int, chosen: List[int]) -> None:
        depth = len(chosen)
        needed = n - depth
        for i in range(start, len(pool) - needed + 1):
            current = limit()
            if _worse_than(single_complexity[i], current):
                break  # sorted: no later candidate can do better either
            top_similarity = builder.push(depth, pool_rows[i])
            if current != INF and _PAIR_SAFE > top_similarity > _pair_limit(current):
                continue  # one pair is already too similar
            c = builder.complexity(depth)
            if np.isfinite(c) and not _worse_than(c, current):
                chosen.append(i)
                if needed == 1:
                    found.append((c, tuple(pool[j] for j in chosen)))
                    found.sort()
                    del found[top_n:]
                else:
                    extend(i + 1, chosen)
                chosen.pop()

    extend(0, [])

    suggestions = []
    for c, extra in found:
        matrix = _similarity_submatrix(panel + list(extra), base, mapped_index)
        suggestions.append({
            "fluors": list(extra),
            "matrix_complexity": c,
            "values": [[round(float(v), 4) for v in row] for row in matrix],
        })
    return {
        "base_complexity": base_complexity,
        "candidates": len(pool),
        "combinations": math.comb(len(pool), n),
        "suggestions": suggestions,
    }


def find_best_panels(
    antibodies: Optional[Dict[str, List[str]]] = None,
    csv_path: Optional[str] = None,
    sim_df: Optional[SimilarityMatrix] = None,
    sim_names: Optional[Sequence[str]] = None,
    det_map: Optional[Dict[str, List[str]]] = None,
    sim_path: Optional[str] = None,
    det_path: Optional[str] = None,
    top_n: int = 10,
    backbone: Optional[Dict[str, str]] = None,
    aliases: Optional[Sequence[Tuple[str, str]]] = None,
    alias_path: Optional[str] = None,
    max_assigned: Optional[int] = None,
) -> dict:
    """
    Find the exact global top_n panels among all maximum-size feasible panels.
    With `max_assigned`, "maximum-size" means at most that many antibodies.

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
        antibodies, sim_df, sim_names, det_map, top_n, backbone, aliases,
        max_assigned,
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
        top_n=10,
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
