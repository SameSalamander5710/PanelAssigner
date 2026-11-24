"""A panel is ranked by its fluorophores, but antibodies that share fluorophores
can usually be given them in more than one way, all with the same score. These
tests cover `panel_assignments` (every such way) and what `worker.finalize`
sends the page so that it can step through them."""

import itertools
import json
import os
import random
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import Panel_optimized as P


def _load_worker():
    """worker.py reads ./constants at import time, so import it from the repo root."""
    cwd = os.getcwd()
    os.chdir(ROOT)
    try:
        import worker
    finally:
        os.chdir(cwd)
    return worker


def _panels(worker, antibodies, top_n=5, backbone=None, max_free=None):
    """The ranked panels (search_subtree's dicts) and the search data."""
    data = P.prepare_search(
        antibodies, worker.REF_SIM, worker.REF_NAMES, worker.REF_DET, top_n,
        backbone, worker.ALIASES, max_free=max_free, canonical=True,
    )
    found = []
    for task in P.make_search_tasks(data):
        best, _, _ = P.search_subtree(data, task, P.INF)
        found.extend(best)
    found.sort(key=lambda item: (item[0], item[1]))
    return data, [item[2] for item in found[:top_n]]


def _free(panel, assignment):
    """An assignment without the fixed backbone, comparable as a set of pairs."""
    return frozenset((a, f) for a, f in assignment.items() if a not in panel["backbone"])


def _brute_force(data, panel):
    """Every way to give each fluorophore of the panel its own antibody, by trying them all."""
    fluors = sorted(f for a, f in panel["assignments"].items() if a not in panel["backbone"])
    found = set()
    for chosen in itertools.permutations(data["free_original_order"], len(fluors)):
        if all(f in data["options"][a] for a, f in zip(chosen, fluors)):
            found.add(frozenset(zip(chosen, fluors)))
    return found


class PanelAssignmentsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.worker = _load_worker()

    def test_antibodies_sharing_two_fluorophores_can_swap(self):
        data, panels = _panels(self.worker, {
            "CD45": ["FITC", "BV421"], "CD11b": ["FITC", "BV421"],
        })
        self.assertEqual(len(panels), 1)  # one set of fluorophores...
        assignments, more = P.panel_assignments(data, panels[0])
        self.assertFalse(more)
        self.assertEqual(len(assignments), 2)  # ...dealt out in two ways
        self.assertEqual(
            {tuple(sorted(a.items())) for a in assignments},
            {(("CD11b", "BV421"), ("CD45", "FITC")), (("CD11b", "FITC"), ("CD45", "BV421"))},
        )
        self.assertEqual(assignments[0], panels[0]["assignments"])  # what the table always showed

    def test_one_assignment_when_there_is_no_choice(self):
        data, panels = _panels(self.worker, {"CD3": ["FITC"], "CD4": ["BV421"]})
        assignments, more = P.panel_assignments(data, panels[0])
        self.assertEqual(assignments, [panels[0]["assignments"]])
        self.assertFalse(more)

    def test_matches_trying_every_assignment(self):
        rng = random.Random(11)
        pool = list(self.worker.REF_NAMES)[:60]
        panels_checked = with_choice = 0
        for _ in range(150):
            n_fluors = rng.randint(2, 6)
            fluors = rng.sample(pool, n_fluors)
            antibodies = {
                f"Ab{i}": rng.sample(fluors, rng.randint(1, min(n_fluors, 4)))
                for i in range(rng.randint(2, 5))
            }
            backbone = {"Fixed": rng.choice([f for f in pool if f not in fluors])} if rng.random() < 0.3 else None
            data, panels = _panels(
                self.worker, antibodies, top_n=4, backbone=backbone,
                max_free=rng.choice([None, None, 1, 2, 3]),
            )
            for panel in panels:
                assignments, more = P.panel_assignments(data, panel, limit=100_000)
                got = [_free(panel, a) for a in assignments]
                self.assertFalse(more)
                self.assertEqual(len(got), len(set(got)), "an assignment is listed twice")
                self.assertEqual(set(got), _brute_force(data, panel), antibodies)
                self.assertEqual(assignments[0], panel["assignments"])
                for assignment in assignments:
                    self.assertEqual(  # the fixed backbone is in every one, unchanged
                        {a: f for a, f in assignment.items() if a in panel["backbone"]}, panel["backbone"])
                    self.assertEqual(  # keys in input order, like panel["assignments"]
                        list(assignment), [a for a in data["original_order"] if a in assignment])
                panels_checked += 1
                with_choice += len(assignments) > 1
        self.assertGreater(panels_checked, 100)
        self.assertGreater(with_choice, 30)  # the comparison is not only about single answers

    def test_order_follows_the_input(self):
        # three antibodies that can take any of three fluorophores, listed in this order
        fluors = ["FITC", "BV421", "APC"]
        data, panels = _panels(self.worker, {f"Ab{i}": list(fluors) for i in range(3)}, top_n=1)
        assignments, more = P.panel_assignments(data, panels[0])
        self.assertFalse(more)
        self.assertEqual(len(assignments), 6)
        base = panels[0]["assignments"]
        rest = [tuple(a[f"Ab{i}"] for i in range(3)) for a in assignments[1:]]
        # the others: Ab0's first-listed fluorophore first, then Ab1's, ...
        expected = [p for p in itertools.permutations(fluors) if p != tuple(base[f"Ab{i}"] for i in range(3))]
        expected.sort(key=lambda p: [fluors.index(f) for f in p])
        self.assertEqual(rest, expected)

    def test_limit_and_more(self):
        fluors = list(self.worker.REF_NAMES)[:6]
        data, panels = _panels(self.worker, {f"Ab{i}": list(fluors) for i in range(6)}, top_n=1)
        panel = panels[0]
        everything, more = P.panel_assignments(data, panel, limit=720)  # 6! ways
        self.assertEqual((len(everything), more), (720, False))
        self.assertEqual(len({tuple(a.items()) for a in everything}), 720)

        some, more = P.panel_assignments(data, panel, limit=50)
        self.assertEqual((len(some), more), (50, True))
        self.assertEqual(some[0], panel["assignments"])
        self.assertEqual(len({tuple(a.items()) for a in some}), 50)

        one, more = P.panel_assignments(data, panel, limit=1)
        self.assertEqual((one, more), ([panel["assignments"]], True))

        default, more = P.panel_assignments(data, panel)
        self.assertEqual((len(default), more), (P.MAX_ALTERNATIVES, True))
        with self.assertRaises(ValueError):
            P.panel_assignments(data, panel, limit=0)

    def test_limited_panels_also_vary_which_antibodies_are_left_out(self):
        data, panels = _panels(self.worker, {
            "A": ["FITC", "BV421"], "B": ["FITC", "BV421"], "C": ["FITC", "BV421"],
        }, max_free=2)
        self.assertEqual(len(panels), 1)
        assignments, _ = P.panel_assignments(data, panels[0])
        self.assertEqual(len(assignments), 6)  # 3 antibodies, 2 places
        self.assertEqual({tuple(sorted(a)) for a in assignments}, {("A", "B"), ("A", "C"), ("B", "C")})
        self.assertTrue(all(len(a) == 2 for a in assignments))

    def test_backbone_is_the_same_in_every_assignment(self):
        data, panels = _panels(self.worker, {
            "CD45": ["FITC", "BV421"], "CD11b": ["FITC", "BV421"],
        }, backbone={"Viability": "APC"})
        assignments, _ = P.panel_assignments(data, panels[0])
        self.assertEqual(len(assignments), 2)
        for assignment in assignments:
            self.assertEqual(list(assignment)[0], "Viability")
            self.assertEqual(assignment["Viability"], "APC")

    def test_backbone_only_has_one_assignment(self):
        data = P.prepare_search(
            {}, self.worker.REF_SIM, self.worker.REF_NAMES, self.worker.REF_DET, 1,
            {"Viability": "APC"}, self.worker.ALIASES, canonical=True,
        )
        best, _, _ = P.search_subtree(data, P.make_search_tasks(data)[0], P.INF)
        assignments, more = P.panel_assignments(data, best[0][2])
        self.assertEqual(assignments, [{"Viability": "APC"}])
        self.assertFalse(more)


class FinalizeAssignmentsTest(unittest.TestCase):
    """What the page gets: the stepper under a panel's name and the data to step with."""

    @classmethod
    def setUpClass(cls):
        cls.worker = _load_worker()

    def _compute(self, text, backbone="", top_n=10, max_fit=""):
        info = json.loads(self.worker.setup(text, str(top_n), backbone, max_fit, "", P.DEFAULT_METRIC))
        self.assertTrue(info["ok"], info)
        results = [self.worker.search(i).partition("\n")[2] for i in range(info["n_tasks"])]
        out = json.loads(self.worker.finalize("[" + ",".join(results) + "]"))
        return out, json.loads(out["matrices_json"])

    def _table_column(self, html, panel_no):
        """The cells (class, text) of one panel's column of the results table."""
        rows = re.findall(r"<tr><td class=\"num\">\d+</td><td>.*?</td>(.*?)</tr>", html)
        cells = [re.findall(r'<td class="([^"]*)"[^>]*>(.*?)</td>', row)[panel_no - 1] for row in rows]
        return [(cls, text) for cls, text in cells]

    def test_the_example_from_the_bug_report(self):
        out, shown = self._compute("CD45,FITC\nCD45,BV421\nCD11b,FITC\nCD11b,BV421")
        self.assertEqual(shown["antibodies"], ["CD45", "CD11b"])
        panel = shown["panels"][0]
        self.assertEqual(len(panel["assignments"]), 2)
        self.assertFalse(panel["more_assignments"])
        # each way gives, per fluorophore of the panel, the table row of its antibody
        by_name = {}
        for assignment in panel["assignments"]:
            by_name[tuple(sorted((shown["antibodies"][a], f["name"])
                                 for a, f in zip(assignment, panel["fluors"])))] = True
        self.assertEqual(set(by_name), {
            (("CD11b", "BV421"), ("CD45", "FITC")), (("CD11b", "FITC"), ("CD45", "BV421")),
        })
        # the table shows the first, and its header says there are two
        self.assertEqual(self._table_column(out["html"], 1), [("assigned", "FITC"), ("assigned", "BV421")])
        self.assertIn('<span class="alt-count" aria-live="polite">1/2</span>', out["html"])
        self.assertIn("nocopy", out["html"])
        self.assertEqual(panel["fluors"][0]["antibody"], "CD45")  # unchanged for older readers

    def test_no_stepper_for_a_panel_without_a_choice(self):
        out, shown = self._compute("CD3,FITC\nCD4,BV421\nCD8,APC", top_n=3)
        self.assertNotIn("alt-nav", out["html"])
        self.assertNotIn("alt-btn", out["html"])
        for panel in shown["panels"]:
            self.assertEqual(len(panel["assignments"]), 1)

    def test_stepper_only_under_panels_that_have_a_choice(self):
        out, shown = self._compute("A,FITC\nA,BV421\nB,FITC\nB,BV421\nB,APC", top_n=3)
        counts = [len(p["assignments"]) for p in shown["panels"]]
        self.assertEqual(sorted(counts), [1, 1, 2])
        headers = re.findall(r"<th>(Panel \d+.*?)</th>", out["html"])
        self.assertEqual([("alt-nav" in h) for h in headers], [c > 1 for c in counts])

    def test_more_is_marked_and_counted_as_a_minimum(self):
        fluors = list(self.worker.REF_NAMES)[:6]
        text = "\n".join(f"Ab{i},{f}" for i in range(6) for f in fluors)
        out, shown = self._compute(text, top_n=1)
        panel = shown["panels"][0]
        self.assertEqual(len(panel["assignments"]), P.MAX_ALTERNATIVES)
        self.assertTrue(panel["more_assignments"])
        self.assertIn(f'1/{P.MAX_ALTERNATIVES}+</span>', out["html"])
        self.assertIn(f'data-widest="{P.MAX_ALTERNATIVES}/{P.MAX_ALTERNATIVES}+"', out["html"])
        self.assertIn(f"at least {P.MAX_ALTERNATIVES} ways", out["html"])

    def test_every_assignment_is_valid_and_keeps_the_matrix_order(self):
        text = "A,FITC\nA,BV421\nA,APC\nB,FITC\nB,BV421\nB,APC\nC,FITC\nC,BV421"
        out, shown = self._compute(text, backbone="Viability,PE", top_n=3)
        options = {"A": {"FITC", "BV421", "APC"}, "B": {"FITC", "BV421", "APC"}, "C": {"FITC", "BV421"}}
        for panel in shown["panels"]:
            names = [f["name"] for f in panel["fluors"]]
            seen = set()
            for assignment in panel["assignments"]:
                rows = [shown["antibodies"][a] for a in assignment]
                self.assertEqual(len(set(rows)), len(rows), "an antibody has two fluorophores")
                for antibody, fluor in zip(rows, names):
                    if antibody != "Viability":
                        self.assertIn(fluor, options[antibody])
                self.assertEqual(rows[names.index("PE")], "Viability")  # the backbone stays
                seen.add(tuple(rows))
            self.assertEqual(len(seen), len(panel["assignments"]), "an assignment is listed twice")
            # all of them are for the same fluorophores in the same matrix order
            self.assertEqual(len(panel["values"]), len(names))

    def test_stepper_is_clear_of_html_in_names(self):
        out, _ = self._compute("<b>X</b>,FITC\n<b>X</b>,BV421\nY&Z,FITC\nY&Z,BV421")
        self.assertNotIn("<b>X</b>", out["html"])
        self.assertIn("&lt;b&gt;X&lt;/b&gt;", out["html"])


if __name__ == "__main__":
    unittest.main()
