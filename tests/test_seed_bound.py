"""seed_bound builds its panels incrementally (one matrix row and one augmenting
path per added fluorophore). It must give the same bound as the plain version
that rebuilds the matrix and the matching from names for every candidate."""

import os
import random
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


def reference_seed_bound(data):
    """The straightforward version: every candidate is checked and scored from scratch."""
    afb, fl = data["antibodies_by_fluor"], data["all_fluors"]
    k, top_n = data["free_max_count"], data["top_n"]

    def feasible(items):
        return len(P._maximum_matching(items, afb)) == len(items)

    chosen = []
    while len(chosen) < k:
        pick, pick_c = None, P.INF
        for f in fl:
            if f in chosen or not feasible(chosen + [f]):
                continue
            c = P._set_complexity(chosen + [f], data)
            if pick is None or c < pick_c:
                pick, pick_c = f, c
        chosen.append(pick)

    found = {frozenset(chosen): P._set_complexity(chosen, data)}
    for pos in range(k):
        rest = chosen[:pos] + chosen[pos + 1:]
        for f in fl:
            if f in chosen:
                continue
            cand = rest + [f]
            key = frozenset(cand)
            if key not in found and feasible(cand):
                found[key] = P._set_complexity(cand, data)
    values = sorted(found.values())
    return values[top_n - 1] if len(values) >= top_n else P.INF


class SeedBoundTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.worker = _load_worker()
        cls.names = [
            n for n in cls.worker.REF_NAMES
            if cls.worker.REF_TYPES.get(P.normalize_fluorophore(n), "antibody") == "antibody"
        ]

    def test_matches_the_straightforward_version(self):
        w, rng, checked = self.worker, random.Random(7), 0
        for _ in range(60):
            n_antibodies, n_options = rng.randint(3, 12), rng.randint(1, 6)
            antibodies = {f"Ab{a}": rng.sample(self.names, n_options) for a in range(n_antibodies)}
            if rng.random() < 0.3:  # a fluorophore with no similarity data
                antibodies["Ab0"].append("NotInTheMatrix")
            backbone = {}
            if rng.random() < 0.5:
                backbone = {f"Bb{i}": rng.choice(self.names) for i in range(rng.randint(1, 2))}
            metric = rng.choice(list(P.METRICS))
            top_n = rng.choice([1, 3, 10, 50])
            try:
                data = P.prepare_search(
                    antibodies, w.REF_SIM, w.REF_NAMES, w.REF_DET, top_n,
                    backbone=backbone, aliases=w.ALIASES, canonical=True,
                    max_free=rng.choice([None, None, 3, 6]), metric=metric,
                )
            except ValueError:  # e.g. a backbone fluorophore that is also an option
                continue
            if data is None:
                continue
            expected, got = reference_seed_bound(data), P.seed_bound(data)
            if expected == P.INF:
                self.assertEqual(got, P.INF)
            else:
                self.assertAlmostEqual(got, expected, delta=1e-12 * max(1.0, expected),
                                       msg=f"{metric} top_n={top_n}")
            checked += 1
        self.assertGreater(checked, 30)  # the loop must really have compared something

    def test_backbone_only(self):
        w = self.worker
        data = P.prepare_search(
            {"Ab": [self.names[0]]}, w.REF_SIM, w.REF_NAMES, w.REF_DET, 1,
            backbone={"Bb1": self.names[1], "Bb2": self.names[2]}, aliases=w.ALIASES,
            canonical=True, max_free=None,
        )
        data = dict(data, free_max_count=0, all_fluors=[])  # nothing free to add
        self.assertAlmostEqual(P.seed_bound(data), reference_seed_bound(data), places=12)


if __name__ == "__main__":
    unittest.main()
