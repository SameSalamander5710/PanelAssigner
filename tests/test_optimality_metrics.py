"""Tests for the four panel metrics (condition index, A-, D- and E-optimality).

Checks the closed-form scores, the pruning bounds that make the branch and
bound search exact, and that the optimizer returns exactly what exhaustive
enumeration finds, for every metric.
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
from tests.test_aoptimality_regression import _make_case, _similarity

ALL = list(P.METRICS)


def _random_spd(rng, n):
    """A random positive definite matrix with unit diagonal (a Gram matrix)."""
    spectra = np.abs(rng.normal(size=(n, n + 6))) + 0.05 * rng.normal(size=(n, n + 6))
    spectra /= np.linalg.norm(spectra, axis=1, keepdims=True)
    return spectra @ spectra.T


def _score(matrix, metric):
    values = np.asarray(matrix, dtype=float)
    names = [f"F{i}" for i in range(len(values))]
    return P._complexity_for_set(
        names, values, {name: i for i, name in enumerate(names)}, metric
    )


class ScoreValuesTest(unittest.TestCase):
    def test_two_by_two_closed_forms(self):
        rho = 0.5
        S = [[1.0, rho], [rho, 1.0]]
        self.assertAlmostEqual(_score(S, "condition_index"), math.sqrt((1 + rho) / (1 - rho)), places=12)
        self.assertAlmostEqual(_score(S, "a_optimality"), 1 / (1 - rho * rho), places=12)
        self.assertAlmostEqual(_score(S, "d_optimality"), (1 - rho * rho) ** -0.5, places=12)
        self.assertAlmostEqual(_score(S, "e_optimality"), 1 / (1 - rho), places=12)

    def test_orthogonal_and_single_panels_score_one(self):
        for metric in ALL:
            self.assertEqual(_score(np.eye(5), metric), 1.0, metric)
            self.assertEqual(_score([[1.0]], metric), 1.0, metric)

    def test_singular_and_nonfinite_panels_are_infinite(self):
        for metric in ALL:
            self.assertEqual(_score([[1.0, 1.0], [1.0, 1.0]], metric), math.inf, metric)
            self.assertEqual(_score([[1.0, np.nan], [np.nan, 1.0]], metric), math.inf, metric)

    def test_metrics_are_ordered_and_at_least_one(self):
        # max >= arithmetic mean >= geometric mean >= 1 of the values 1/lambda
        rng = np.random.default_rng(1)
        for _ in range(200):
            S = _random_spd(rng, int(rng.integers(2, 9)))
            e, a, d = (_score(S, m) for m in ("e_optimality", "a_optimality", "d_optimality"))
            self.assertGreaterEqual(e, a - 1e-9)
            self.assertGreaterEqual(a, d - 1e-9)
            self.assertGreaterEqual(d, 1.0 - 1e-9)
            self.assertGreaterEqual(_score(S, "condition_index"), 1.0 - 1e-9)

    def test_d_optimality_matches_determinant(self):
        rng = np.random.default_rng(2)
        S = _random_spd(rng, 6)
        self.assertAlmostEqual(_score(S, "d_optimality"), np.linalg.det(S) ** (-1 / 6), places=9)
        self.assertAlmostEqual(_score(S, "a_optimality"), np.trace(np.linalg.inv(S)) / 6, places=9)
        self.assertAlmostEqual(_score(S, "e_optimality"), np.linalg.eigvalsh(np.linalg.inv(S))[-1], places=9)

    def test_unknown_metric_is_rejected(self):
        with self.assertRaises(ValueError):
            P.get_metric("f_optimality")
        with self.assertRaises(ValueError):
            P.prepare_search({"A": ["F0"]}, _similarity(np.random.default_rng(0), ["F0"]),
                             ["F0"], None, 1, metric="nope")


class BoundTest(unittest.TestCase):
    """The pruning bounds must never exceed the true final score."""

    def test_subset_lower_bound_never_exceeds_final_score(self):
        rng = np.random.default_rng(3)
        for metric in ALL:
            scorer = P.get_metric(metric)
            for _ in range(300):
                n = int(rng.integers(2, 9))
                S = _random_spd(rng, n)
                m = int(rng.integers(1, n))
                subset = rng.choice(n, size=m, replace=False)
                sub_score = _score(S[np.ix_(subset, subset)], metric)
                bound = scorer.lower_bound(sub_score, m, n) if m >= 2 else 1.0
                self.assertLessEqual(bound, _score(S, metric) * (1 + 1e-9), metric)

    def test_pair_limit_is_sound(self):
        rng = np.random.default_rng(4)
        for metric in ALL:
            scorer = P.get_metric(metric)
            for _ in range(300):
                n = int(rng.integers(2, 9))
                S = _random_spd(rng, n)
                top = np.abs(S - np.eye(n)).max()
                final = _score(S, metric)
                for bound in np.linspace(1.0, final * 1.5, 12):
                    if top > scorer.pair_limit(float(bound), n):
                        self.assertGreater(final, bound * (1 - 1e-9), metric)

    def test_pair_limit_is_not_vacuous(self):
        # a pair this similar must be able to rule out a panel of 10 under a good limit
        for metric in ALL:
            self.assertLess(P.get_metric(metric).pair_limit(2.0, 10), 1.0, metric)


def _exhaustive_panels(antibodies, sim, names, top_n, metric, backbone=None, max_assigned=None):
    data = P.prepare_search(
        antibodies, sim, names, None, top_n,
        backbone=backbone, max_assigned=max_assigned, metric=metric,
    )
    fixed = list(data["backbone_fluors"])
    k = data["free_max_count"]
    found = []
    for combo in itertools.combinations(data["all_fluors"], k):
        if len(P._maximum_matching(list(combo), data["antibodies_by_fluor"])) != k:
            continue
        score = P._complexity_for_set(fixed + list(combo), data["base"], data["mapped_index"], metric)
        found.append((score, tuple(sorted(fixed + list(combo)))))
    found.sort(key=lambda item: (item[0], item[1]))
    return found[:top_n]


class ExhaustiveEquivalenceTest(unittest.TestCase):
    def assert_same(self, expected, got, msg=""):
        self.assertEqual(len(expected), len(got), msg)
        for (e_score, e_set), (g_score, g_set) in zip(sorted(expected), sorted(got)):
            self.assertAlmostEqual(e_score, g_score, delta=1e-9, msg=msg)
            self.assertEqual(e_set, g_set, msg)

    def test_panels_match_exhaustive_for_every_metric(self):
        for offset, metric in enumerate(ALL):
            rng = np.random.default_rng(20261120 + offset)
            for trial in range(60):
                antibodies, sim, names = _make_case(
                    rng, int(rng.integers(2, 6)), int(rng.integers(3, 8)))
                top_n = int(rng.integers(1, 4))
                expected = _exhaustive_panels(antibodies, sim, names, top_n, metric)
                result = P.find_best_panels(
                    antibodies=antibodies, sim_df=sim, sim_names=names,
                    top_n=top_n, metric=metric)
                self.assertEqual(result["metric"], metric)
                got = [(p["matrix_complexity"], tuple(sorted(p["assignments"].values())))
                       for p in result["ranked_panels"]]
                self.assert_same(expected, got, f"{metric} trial {trial}")

    def test_backbone_and_size_limit_match_exhaustive_for_every_metric(self):
        for offset, metric in enumerate(ALL):
            rng = np.random.default_rng(20261130 + offset)
            for trial in range(40):
                antibodies, sim, names = _make_case(
                    rng, int(rng.integers(2, 6)), int(rng.integers(4, 9)))
                first = next(iter(antibodies))
                backbone = {first: antibodies[first][0]}
                max_assigned = int(rng.integers(1, 5))
                top_n = int(rng.integers(1, 4))
                expected = _exhaustive_panels(
                    antibodies, sim, names, top_n, metric, backbone, max_assigned)
                result = P.find_best_panels(
                    antibodies=antibodies, sim_df=sim, sim_names=names, top_n=top_n,
                    backbone=backbone, max_assigned=max_assigned, metric=metric)
                got = [(p["matrix_complexity"], tuple(sorted(p["assignments"].values())))
                       for p in result["ranked_panels"]]
                self.assert_same(expected, got, f"{metric} trial {trial}")

    def test_suggestions_match_exhaustive_for_every_metric(self):
        for offset, metric in enumerate(ALL):
            rng = np.random.default_rng(20261140 + offset)
            for trial in range(60):
                names = [f"F{i}" for i in range(int(rng.integers(4, 9)))]
                sim = _similarity(rng, names)
                panel = [str(f) for f in rng.choice(
                    names, size=int(rng.integers(0, min(4, len(names)))), replace=False)]
                available = [x for x in names if x not in panel]
                n = int(rng.integers(1, min(3, len(available)) + 1))

                base = sim.symmetrized()
                sim_map = {name: i for i, name in enumerate(sim.index)}
                expected = []
                for combo in itertools.combinations(available, n):
                    score = P._complexity_for_set(panel + list(combo), base, sim_map, metric)
                    if np.isfinite(score):
                        expected.append((score, tuple(sorted(combo))))
                expected.sort(key=lambda item: (item[0], item[1]))

                result = P.suggest_additions(sim, names, None, panel, n, top_n=5, metric=metric)
                self.assertEqual(result["metric"], metric)
                got = [(i["matrix_complexity"], tuple(sorted(i["fluors"])))
                       for i in result["suggestions"]]
                self.assert_same(expected[:5], got, f"{metric} trial {trial}")


if __name__ == "__main__":
    unittest.main()
