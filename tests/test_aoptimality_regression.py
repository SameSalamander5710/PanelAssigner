"""Regression tests for the A-optimality scoring and branch-and-bound search.

The optimizer is compared against exhaustive enumeration on small, deterministic
random inputs. The generated similarity matrices are strictly diagonally
dominant, so they are positive definite and every panel has a finite score.
"""

import itertools
import math
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import Panel_optimized as P


def _similarity(rng, names):
    """A random symmetric, strictly diagonally dominant similarity matrix."""
    n = len(names)
    values = rng.uniform(-0.6, 0.9, size=(n, n))
    values = (values + values.T) / 2.0
    np.fill_diagonal(values, 0.0)
    row_abs = np.abs(values).sum(axis=1).max()
    if row_abs >= 0.95:
        values *= 0.95 / row_abs
    return P.SimilarityMatrix(names, names, values + np.eye(n))


def _make_case(rng, n_antibodies, n_fluorophores, max_options=3):
    names = [f"F{i}" for i in range(n_fluorophores)]
    antibodies = {}
    for a in range(n_antibodies):
        size = int(rng.integers(1, min(max_options, n_fluorophores) + 1))
        options = rng.choice(names, size=size, replace=False)
        antibodies[f"A{a}"] = [str(f) for f in options]
    return antibodies, _similarity(rng, names), names


def _exhaustive_panels(
    antibodies, sim, names, top_n, backbone=None, max_assigned=None
):
    data = P.prepare_search(
        antibodies, sim, names, None, top_n,
        backbone=backbone, max_assigned=max_assigned,
    )
    if data is None:
        return [], None

    fixed = list(data["backbone_fluors"])
    k = data["free_max_count"]
    antibodies_by_fluor = data["antibodies_by_fluor"]
    found = []
    for combo in itertools.combinations(data["all_fluors"], k):
        if len(P._maximum_matching(list(combo), antibodies_by_fluor)) != k:
            continue
        score = P._complexity_for_set(
            fixed + list(combo), data["base"], data["mapped_index"]
        )
        found.append((score, tuple(sorted(fixed + list(combo)))))
    found.sort(key=lambda item: (item[0], item[1]))
    return found[:top_n], data


def _exhaustive_suggestions(sim, names, panel, n, top_n=5):
    taken = {P.normalize_fluorophore(name) for name in panel}
    candidates = [
        name for name in dict.fromkeys(names)
        if P.normalize_fluorophore(name) not in taken
    ]
    base = sim.symmetrized()
    sim_map = {name: i for i, name in enumerate(sim.index)}

    found = []
    for combo in itertools.combinations(candidates, n):
        score = P._complexity_for_set(
            list(panel) + list(combo), base, sim_map
        )
        if np.isfinite(score):
            found.append((score, tuple(sorted(combo))))
    found.sort(key=lambda item: (item[0], item[1]))
    return found[:top_n]


class ExhaustiveVersusOptimizedTest(unittest.TestCase):
    def assert_panels_match(self, expected, result):
        got = [
            (panel["matrix_complexity"],
             tuple(sorted(panel["assignments"].values())))
            for panel in result["ranked_panels"]
        ]
        got.sort(key=lambda item: (item[0], item[1]))
        expected = sorted(expected, key=lambda item: (item[0], item[1]))

        self.assertEqual(len(got), len(expected))
        for (expected_score, expected_set), (got_score, got_set) in zip(expected, got):
            self.assertAlmostEqual(
                expected_score, got_score, delta=1e-9,
                msg=f"{expected_set} vs {got_set}",
            )
            self.assertEqual(expected_set, got_set)

    def test_random_panels_match_exhaustive_enumeration(self):
        rng = np.random.default_rng(20261113)
        for trial in range(100):
            antibodies, sim, names = _make_case(
                rng,
                n_antibodies=int(rng.integers(2, 6)),
                n_fluorophores=int(rng.integers(3, 8)),
            )
            top_n = int(rng.integers(1, 4))

            expected, data = _exhaustive_panels(antibodies, sim, names, top_n)
            self.assertIsNotNone(data, f"trial {trial}: nothing to search")
            result = P.find_best_panels(
                antibodies=antibodies, sim_df=sim, sim_names=names, top_n=top_n
            )
            self.assert_panels_match(expected, result)

    def test_random_panels_with_backbone_and_size_limit_match_exhaustive(self):
        rng = np.random.default_rng(20261114)
        for trial in range(50):
            antibodies, sim, names = _make_case(
                rng,
                n_antibodies=int(rng.integers(2, 6)),
                n_fluorophores=int(rng.integers(4, 9)),
            )
            backbone_antibody = next(iter(antibodies))
            backbone = {backbone_antibody: antibodies[backbone_antibody][0]}
            max_assigned = int(rng.integers(1, 5))
            top_n = int(rng.integers(1, 4))

            expected, data = _exhaustive_panels(
                antibodies, sim, names, top_n,
                backbone=backbone, max_assigned=max_assigned,
            )
            self.assertIsNotNone(data, f"trial {trial}: nothing to search")
            result = P.find_best_panels(
                antibodies=antibodies, sim_df=sim, sim_names=names, top_n=top_n,
                backbone=backbone, max_assigned=max_assigned,
            )
            self.assert_panels_match(expected, result)

    def test_random_suggestions_match_exhaustive_enumeration(self):
        rng = np.random.default_rng(20261115)
        for trial in range(100):
            names = [f"F{i}" for i in range(int(rng.integers(4, 9)))]
            sim = _similarity(rng, names)
            panel_size = int(rng.integers(0, min(4, len(names))))
            panel = [str(f) for f in rng.choice(names, size=panel_size, replace=False)]
            available = [name for name in names if name not in panel]
            n = int(rng.integers(1, min(3, len(available)) + 1))

            expected = _exhaustive_suggestions(sim, names, panel, n)
            result = P.suggest_additions(
                sim, names, None, panel, n, top_n=5
            )

            self.assertEqual(result["candidates"], len(available))
            self.assertEqual(result["combinations"], math.comb(len(available), n))
            got = [
                (item["matrix_complexity"], tuple(sorted(item["fluors"])))
                for item in result["suggestions"]
            ]
            got.sort(key=lambda item: (item[0], item[1]))
            expected = sorted(expected, key=lambda item: (item[0], item[1]))

            self.assertEqual(len(got), len(expected), f"trial {trial}")
            for (expected_score, expected_set), (got_score, got_set) in zip(expected, got):
                self.assertAlmostEqual(
                    expected_score, got_score, delta=1e-9,
                    msg=f"trial {trial}: {expected_set} vs {got_set}",
                )
                self.assertEqual(expected_set, got_set)


class AOptimalityMetricTest(unittest.TestCase):
    def score(self, matrix):
        values = np.asarray(matrix, dtype=float)
        names = [f"F{i}" for i in range(len(values))]
        return P._complexity_for_set(
            names, values, {name: i for i, name in enumerate(names)}
        )

    def test_single_fluorophore_scores_one(self):
        self.assertEqual(self.score([[1.0]]), 1.0)

    def test_orthogonal_panel_scores_one(self):
        self.assertEqual(self.score(np.eye(5)), 1.0)

    def test_two_by_two_formula(self):
        rho = 0.5
        self.assertAlmostEqual(
            self.score([[1.0, rho], [rho, 1.0]]),
            1.0 / (1.0 - rho * rho),
            places=12,
        )

    def test_singular_panel_is_infinite(self):
        self.assertEqual(
            self.score([[1.0, 1.0], [1.0, 1.0]]), float("inf")
        )

    def test_nonfinite_panel_is_infinite(self):
        self.assertEqual(
            self.score([[1.0, np.nan], [np.nan, 1.0]]), float("inf")
        )


if __name__ == "__main__":
    unittest.main()
