"""The "Load example panel" files must stay usable.

constants/Example_Panel.csv and constants/Example_Backbone.csv are what the
WebUI puts in the two text boxes. They should read cleanly (no substituted or
unmapped fluorophores), show off aliases, and give a panel that the search
can fully assign.
"""

import json
import os
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import Panel_optimized as P

PANEL = (ROOT / "constants" / "Example_Panel.csv").read_text(encoding="utf-8").strip()
BACKBONE = (ROOT / "constants" / "Example_Backbone.csv").read_text(encoding="utf-8").strip()


def _load_worker():
    """worker.py reads ./constants at import time, so import it from the repo root."""
    cwd = os.getcwd()
    os.chdir(ROOT)
    try:
        import worker
    finally:
        os.chdir(cwd)
    return worker


class ExamplePanelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.worker = _load_worker()
        cls.info = json.loads(cls.worker.summarize(PANEL, BACKBONE, ""))

    def test_files_are_plain_text_without_carriage_returns(self):
        for text in (PANEL, BACKBONE):
            self.assertNotIn("\r", text)
            self.assertTrue(text.startswith("Antibody,Fluorophore"))

    def test_size_and_clean_reading(self):
        info = self.info
        self.assertTrue(info["ok"], info)
        self.assertTrue(10 <= info["antibodies"] <= 14, info["antibodies"])
        self.assertGreaterEqual(info["pairs"], 2 * info["antibodies"])  # real choices
        self.assertGreaterEqual(info["fixed"], 1)  # the viability dye
        self.assertEqual(info["substituted"], [])
        self.assertEqual(info["unmapped"], [])

    def test_row_highlighting_lines_up_with_the_text(self):
        # the page only colours rows when there is exactly one status per line
        self.assertEqual(len(self.info["input_rows"]), len(PANEL.split("\n")))
        self.assertEqual(len(self.info["backbone_rows"]), len(BACKBONE.split("\n")))

    def test_shows_aliases(self):
        worker = self.worker
        known = set(worker.REF_NAMES)
        written = {f for options in worker.parse_pasted(PANEL).values() for f in options}
        written |= set(worker.parse_backbone(BACKBONE).values())
        spelled_differently = [f for f in written if f not in known]
        self.assertGreaterEqual(len(spelled_differently), 5, sorted(spelled_differently))

    def test_every_antibody_gets_a_finite_panel(self):
        worker = self.worker
        antibodies = worker.parse_pasted(PANEL)
        backbone = worker.parse_backbone(BACKBONE)
        result = P.find_best_panels(
            antibodies=antibodies, sim_df=worker.SIM_DF, sim_names=worker.SIM_NAMES,
            det_map=worker.DET_MAP, top_n=3, backbone=backbone, aliases=worker.ALIASES,
        )
        self.assertEqual(result["max_assigned"], result["total_antibodies"])
        self.assertTrue(result["ranked_panels"])
        for panel in result["ranked_panels"]:
            self.assertLess(panel["matrix_complexity"], float("inf"))
            self.assertEqual(len(panel["assignments"]), result["total_antibodies"])


if __name__ == "__main__":
    unittest.main()
