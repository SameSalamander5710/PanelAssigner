import json
import os
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _load_worker():
    cwd = os.getcwd()
    os.chdir(ROOT)
    try:
        import worker
    finally:
        os.chdir(cwd)
    return worker


class PastedInputTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.worker = _load_worker()

    def test_mixed_commas_and_tabs_are_parsed_together(self):
        text = (
            "Antibody,Fluorophore\n"
            "CD3,BV421\n"
            "CD4\tFITC\n"
            "CD8,APC\n"
            "CD4\tPE\n"
        )
        self.assertEqual(self.worker.parse_pasted(text), {
            "CD3": ["BV421"],
            "CD4": ["FITC", "PE"],
            "CD8": ["APC"],
        })
        self.assertEqual(self.worker._row_fluors(text), [
            (1, "BV421"), (2, "FITC"), (3, "APC"), (4, "PE"),
        ])

    def test_mixed_lines_with_windows_line_endings(self):
        text = "CD3,BV421\r\nCD4\tFITC\r\nCD8,APC\r\n"
        self.assertEqual(self.worker.parse_pasted(text), {
            "CD3": ["BV421"], "CD4": ["FITC"], "CD8": ["APC"],
        })
        self.assertEqual(self.worker._row_fluors(text), [
            (0, "BV421"), (1, "FITC"), (2, "APC"),
        ])

    def test_backbone_mixed_delimiters(self):
        text = (
            "Antibody,Fluorophore\n"
            "Viability,BV421\n"
            "Comp\tFITC\n"
            "PE\n"
        )
        self.assertEqual(self.worker.parse_backbone(text), {
            "Viability": "BV421", "Comp": "FITC", "Avoid_1": "PE",
        })
        self.assertEqual(self.worker._row_fluors(text, True), [
            (1, "BV421"), (2, "FITC"), (3, "PE"),
        ])

    def test_quoted_separators_stay_inside_the_value(self):
        text = 'CD3,"BV421,FITC"\nCD4\t"APC,PE"\nCD8,"BV605\tPE"\n'
        self.assertEqual(self.worker.parse_pasted(text), {
            "CD3": ["BV421,FITC"],
            "CD4": ["APC,PE"],
            "CD8": ["BV605\tPE"],
        })

    def test_summary_and_setup_see_every_mixed_row(self):
        text = "CD3,BV421\nCD4\tFITC\nCD8,APC"
        info = json.loads(self.worker.summarize(text, "", ""))
        self.assertTrue(info["ok"], info)
        self.assertEqual((info["antibodies"], info["pairs"]), (3, 3))
        self.assertEqual(len(info["input_rows"]), 3)

        setup = json.loads(self.worker.setup(text, "3", "", "", ""))
        self.assertTrue(setup["ok"], setup)
        self.assertEqual(self.worker.STATE["free_total"], 3)


if __name__ == "__main__":
    unittest.main()
