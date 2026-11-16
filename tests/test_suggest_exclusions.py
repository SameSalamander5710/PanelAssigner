"""Leaving fluorophores out of the suggestions (the "leave out" box of the WebUI).

An exclusion removes fluorophores from the pool of candidates; the best sets are
then the best of what is left, exactly as if those fluorophores did not exist.
Names match ignoring case and spaces, `*` is a wildcard, and every spelling
that an alias rule gives a fluorophore counts ("BV421" is "Brilliant Violet 421").
"""

import os
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import Panel_optimized as P

ALIASES = P.load_fluorophore_aliases(str(ROOT / "constants" / "Fluorophore_Aliases.csv"))
NAMES = ["BV421", "BV605", "BUV395", "PE", "PE-Cy7", "FITC", "APC", "Alexa Fluor 700", "Qdot 605", "PerCP"]


def similarity(names, seed=7):
    """A random symmetric, strictly diagonally dominant matrix (so every set has a finite score)."""
    rng = np.random.default_rng(seed)
    n = len(names)
    values = rng.uniform(0.0, 0.9, size=(n, n))
    values = (values + values.T) / 2.0
    np.fill_diagonal(values, 0.0)
    values *= 0.9 / np.abs(values).sum(axis=1).max()
    return P.SimilarityMatrix(names, names, values + np.eye(n))


def without(sim, names, removed):
    keep = [i for i, name in enumerate(names) if name not in removed]
    kept = [names[i] for i in keep]
    return P.SimilarityMatrix(kept, kept, sim.values[np.ix_(keep, keep)]), kept


def suggest(exclude, panel=("APC",), n=2, names=NAMES, sim=None):
    sim = sim or similarity(names)
    return P.suggest_additions(sim, names, None, list(panel), n, top_n=5, aliases=ALIASES, exclude=exclude)


class ParseExclusionsTest(unittest.TestCase):
    def test_separators_and_repeats(self):
        self.assertEqual(
            P.parse_exclusions("BV*, PE-Cy7\nbv*;  Qdot 605\t\n\n"),
            ["BV*", "PE-Cy7", "Qdot 605"],
        )

    def test_repeat_ignores_case_and_spacing(self):
        self.assertEqual(P.parse_exclusions("Alexa Fluor 700\nalexafluor700"), ["Alexa Fluor 700"])

    def test_empty(self):
        self.assertEqual(P.parse_exclusions(""), [])
        self.assertEqual(P.parse_exclusions(None), [])
        self.assertEqual(P.parse_exclusions(" , ;\n"), [])


class MatchingTest(unittest.TestCase):
    def excluded(self, pattern):
        kept, info = P._apply_exclusions(NAMES, [pattern], ALIASES)
        self.assertEqual(info[0]["removed"], len(NAMES) - len(kept))
        return [name for name in NAMES if name not in kept]

    def test_prefix_wildcard(self):
        self.assertEqual(self.excluded("BV*"), ["BV421", "BV605"])

    def test_wildcard_reaches_alias_spellings(self):
        self.assertEqual(self.excluded("Brilliant Violet*"), ["BV421", "BV605"])
        self.assertEqual(self.excluded("PE-Cyanine*"), ["PE-Cy7"])

    def test_exact_name_in_any_spelling(self):
        self.assertEqual(self.excluded("bv 421"), ["BV421"])
        self.assertEqual(self.excluded("Brilliant Violet 421"), ["BV421"])
        self.assertEqual(self.excluded("AF700"), ["Alexa Fluor 700"])
        self.assertEqual(self.excluded("PE-Cyanine 7"), ["PE-Cy7"])

    def test_wildcard_anywhere(self):
        self.assertEqual(self.excluded("*605"), ["BV605", "Qdot 605"])
        self.assertEqual(self.excluded("PE*"), ["PE", "PE-Cy7", "PerCP"])  # PerCP starts with Pe too

    def test_plain_name_is_not_a_prefix(self):
        self.assertEqual(self.excluded("PE"), ["PE"])  # not PE-Cy7
        self.assertEqual(self.excluded("BV"), [])

    def test_unmatched_pattern_reports_zero(self):
        kept, info = P._apply_exclusions(NAMES, ["Nonsense", "BV*"], ALIASES)
        self.assertEqual([(e["pattern"], e["removed"], e["names"]) for e in info],
                         [("Nonsense", 0, []), ("BV*", 2, ["BV421", "BV605"])])
        self.assertEqual(len(kept), len(NAMES) - 2)


class SuggestWithExclusionsTest(unittest.TestCase):
    def test_nothing_excluded_changes_nothing(self):
        plain = P.suggest_additions(similarity(NAMES), NAMES, None, ["APC"], 2, top_n=5, aliases=ALIASES)
        same = suggest([])
        self.assertEqual(plain["suggestions"], same["suggestions"])
        self.assertEqual(same["excluded"], [])

    def test_equals_a_reference_without_those_fluorophores(self):
        sim = similarity(NAMES)
        for pattern, gone in (("BV*", {"BV421", "BV605"}), ("Qdot 605", {"Qdot 605"}), ("PE*", {"PE", "PE-Cy7", "PerCP"})):
            reduced, kept = without(sim, NAMES, gone)
            expected = P.suggest_additions(reduced, kept, None, ["APC"], 2, top_n=5, aliases=ALIASES)
            got = suggest([pattern], sim=sim)
            self.assertEqual(
                [(s["fluors"], round(s["matrix_complexity"], 9)) for s in got["suggestions"]],
                [(s["fluors"], round(s["matrix_complexity"], 9)) for s in expected["suggestions"]],
                pattern,
            )
            self.assertEqual(got["candidates"], expected["candidates"])
            self.assertEqual(got["combinations"], expected["combinations"])

    def test_nothing_excluded_is_suggested(self):
        got = suggest(["BV*", "Qdot 605"], n=3)
        suggested = {f for s in got["suggestions"] for f in s["fluors"]}
        self.assertFalse(suggested & {"BV421", "BV605", "Qdot 605"})
        self.assertEqual([(e["pattern"], e["names"]) for e in got["excluded"]],
                         [("BV*", ["BV421", "BV605"]), ("Qdot 605", ["Qdot 605"])])

    def test_panel_fluorophores_are_never_filtered(self):
        sim = similarity(NAMES)
        panel = ["BV421", "APC"]
        plain = P.suggest_additions(sim, NAMES, None, panel, 1, top_n=5, aliases=ALIASES)
        got = P.suggest_additions(sim, NAMES, None, panel, 1, top_n=5, aliases=ALIASES, exclude=["BV*"])
        self.assertEqual(plain["base_complexity"], got["base_complexity"])
        self.assertEqual([(e["pattern"], e["names"]) for e in got["excluded"]],
                         [("BV*", ["BV605"])])  # BV421 is in the panel, so it was never a candidate
        self.assertTrue(all(s["fluors"][0] != "BV605" for s in got["suggestions"]))

    def test_too_few_left(self):
        with self.assertRaisesRegex(ValueError, "after leaving out"):
            suggest(["*"], n=1)
        with self.assertRaisesRegex(ValueError, "Only 0"):
            suggest(["*"], n=2)

    def test_pool_count_follows_the_exclusions(self):
        sim = similarity(NAMES)
        args = (sim, NAMES, None, ["APC"])
        all_of_them = P.eligible_fluorophore_count(*args, aliases=ALIASES)
        fewer = P.eligible_fluorophore_count(*args, aliases=ALIASES, exclude=["BV*"])
        self.assertEqual(all_of_them - fewer, 2)


class WorkerSuggestTest(unittest.TestCase):
    """worker.suggest: the text of the box goes in, what was applied comes back."""

    @classmethod
    def setUpClass(cls):
        cwd = os.getcwd()
        os.chdir(ROOT)  # worker.py reads ./constants at import time
        try:
            import worker
        finally:
            os.chdir(cwd)
        cls.worker = worker
        cls.panel = ["BUV395", "BV650", "BUV805", "PE-Cy7", "FITC", "PE", "APC-Fire 750"]

    def setUp(self):
        self.worker.STATE["panel_labels"] = [self.panel]
        self.worker.STATE["metric"] = "a_optimality"

    def run_suggest(self, exclude):
        import json
        return json.loads(self.worker.suggest(1, 2, exclude))

    def test_no_exclusions(self):
        result = self.run_suggest("")
        self.assertTrue(result["ok"], result)
        self.assertEqual(result["exclude"], [])
        self.assertEqual(result["excluded"], [])

    def test_wildcard_with_aliases(self):
        free = self.run_suggest("")
        first = free["suggestions"][0]["fluors"][0]
        prefix = first[:2]  # e.g. "RB" of "RB705": a wildcard that takes the best pick away
        result = self.run_suggest(f"{prefix}*\nBrilliant Violet*, nonsense")
        self.assertTrue(result["ok"], result)
        picked = [f for s in result["suggestions"] for f in s["fluors"]]
        self.assertFalse([f for f in picked if f.lower().startswith(prefix.lower()) or f.startswith("BV")], picked)
        self.assertEqual(result["exclude"], [f"{prefix}*", "Brilliant Violet*", "nonsense"])
        removed = {e["pattern"]: e for e in result["excluded"]}
        self.assertIn(first, removed[f"{prefix}*"]["names"])
        self.assertGreater(removed["Brilliant Violet*"]["removed"], 0)  # an alias spelling
        self.assertTrue(all(n.startswith("BV") for n in removed["Brilliant Violet*"]["names"]))
        self.assertEqual(removed["nonsense"]["removed"], 0)
        self.assertLess(result["candidates"], free["candidates"])

    def test_error_when_nothing_is_left(self):
        result = self.run_suggest("*")
        self.assertFalse(result["ok"])
        self.assertIn("left", result["error"])


if __name__ == "__main__":
    unittest.main()
