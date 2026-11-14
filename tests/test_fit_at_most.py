"""\"Fit at most N\" counts the antibodies of the main box, not the fixed ones.

In the engine that is `max_free`; `max_assigned` (panel size, backbone
included) is unchanged. The WebUI uses `max_free`.
"""

import json
import os
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import Panel_optimized as P
from tests.test_aoptimality_regression import _make_case


def _ranked(result):
    return [(p["matrix_complexity"], tuple(sorted(p["assignments"].items())))
            for p in result["ranked_panels"]]


class MaxFreeTest(unittest.TestCase):
    def test_max_free_equals_max_assigned_with_backbone_added(self):
        for metric in P.METRICS:
            rng = np.random.default_rng(20261201)
            for trial in range(40):
                antibodies, sim, names = _make_case(
                    rng, int(rng.integers(3, 6)), int(rng.integers(4, 9)))
                fixed = list(antibodies)[: int(rng.integers(0, 3))]
                backbone = {a: antibodies[a][0] for a in fixed}
                k = int(rng.integers(1, 5))
                common = dict(antibodies=antibodies, sim_df=sim, sim_names=names,
                              top_n=3, backbone=backbone, metric=metric)
                by_free = P.find_best_panels(max_free=k, **common)
                by_size = P.find_best_panels(max_assigned=len(backbone) + k, **common)
                self.assertEqual(_ranked(by_free), _ranked(by_size), f"{metric} trial {trial}")
                self.assertEqual(by_free["max_assigned"], by_size["max_assigned"])

    def test_the_fixed_backbone_never_uses_up_the_limit(self):
        rng = np.random.default_rng(7)
        antibodies, sim, names = _make_case(rng, 6, 10)
        first = next(iter(antibodies))
        backbone = {first: antibodies[first][0]}
        result = P.find_best_panels(antibodies=antibodies, sim_df=sim, sim_names=names,
                                    top_n=2, backbone=backbone, max_free=2)
        # 2 searched antibodies plus the fixed one (at most: some may not fit)
        self.assertLessEqual(result["max_assigned"], 3)
        for panel in result["ranked_panels"]:
            self.assertIn(first, panel["assignments"])

    def test_both_limits_or_a_limit_below_one_are_rejected(self):
        rng = np.random.default_rng(8)
        antibodies, sim, names = _make_case(rng, 3, 5)
        common = dict(antibodies=antibodies, sim_df=sim, sim_names=names, top_n=1)
        with self.assertRaises(ValueError):
            P.find_best_panels(max_free=2, max_assigned=3, **common)
        with self.assertRaises(ValueError):
            P.find_best_panels(max_free=0, **common)


class WorkerFitAtMostTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cwd = os.getcwd()
        os.chdir(ROOT)  # worker.py reads ./constants at import time
        try:
            import worker
        finally:
            os.chdir(cwd)
        cls.worker = worker
        cls.panel = (ROOT / "constants" / "Example_Panel.csv").read_text(encoding="utf-8")
        cls.backbone = (ROOT / "constants" / "Example_Backbone.csv").read_text(encoding="utf-8")

    def _setup(self, max_fit):
        info = json.loads(self.worker.setup(self.panel, "3", self.backbone, max_fit, ""))
        self.assertTrue(info["ok"], info)
        return info

    def _compute(self, max_fit):
        info = self._setup(max_fit)
        bound, results = float(info["bound"]), []
        for i in range(info["n_tasks"]):
            results.append(self.worker.search(i, repr(bound)).partition("\n")[2])
        return json.loads(self.worker.finalize("[" + ",".join(results) + "]"))

    def test_limit_counts_only_the_main_box(self):
        # the example has 12 antibodies in the main box and one fixed viability dye
        self._setup("")
        self.assertEqual(self.worker.STATE["free_total"], 12)
        self.assertEqual(self.worker.STATE["total"], 13)  # the table still lists the dye
        self.assertEqual(self.worker.STATE["max_count"], 13)

        self._setup("12")  # "12 of 12" is the same as no limit
        self.assertEqual(self.worker.STATE["max_count"], 13)

        res = self._compute("9")
        self.assertEqual(self.worker.STATE["max_count"], 10)  # 9 searched + 1 fixed
        shown = json.loads(res["matrices_json"])
        self.assertEqual(shown["free_total"], 12)
        for panel in shown["panels"]:
            searched = [f for f in panel["fluors"] if not f.get("backbone")]
            self.assertEqual(len(searched), 9)
            self.assertEqual(len(panel["fluors"]), 10)

    def test_zero_is_reported(self):
        info = json.loads(self.worker.setup(self.panel, "3", self.backbone, "0", ""))
        self.assertFalse(info["ok"])
        self.assertIn("at least 1", info["error"])


if __name__ == "__main__":
    unittest.main()
