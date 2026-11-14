"""The live check of the inputs (`refresh_live` in main.py) must never stay stale.

Typing again while a check is still with a worker used to start a second call
on the same worker, which was never answered: the first result was dropped as
outdated and the summary stayed greyed out. These tests load the real main.py
with a stand-in for the `pyscript` module and a worker that, like the browser,
answers one call at a time and never answers a call made while another runs.
"""

import asyncio
import json
import sys
import types
import unittest
from pathlib import Path

MAIN = Path(__file__).resolve().parents[1] / "main.py"
DEBOUNCE = 0.05  # replaces main.LIVE_DEBOUNCE
CALL = 0.5  # how long a worker call takes (longer than the pause before a check)


class _Classes:
    def add(self, *names): pass
    def remove(self, *names): pass


class _Element:
    def __init__(self):
        self.textContent = self.innerText = self.innerHTML = self.value = ""
        self.checked = False
        self.disabled = True
        self.classList = _Classes()

    def getAttribute(self, name):
        return ""


class _Document:
    def __init__(self):
        self.elements = {}
        self.getElementById("topn").value = "10"

    def getElementById(self, name):
        return self.elements.setdefault(name, _Element())


class _Worker:
    def __init__(self):
        self.active = False
        self.overlaps = 0  # calls made while another call was running
        self.fail_next = False
        self.hang_next = False

    async def _call(self, answer):
        if self.active:
            self.overlaps += 1
            await asyncio.Event().wait()  # never answered
        if self.hang_next:
            self.hang_next = False
            await asyncio.Event().wait()
        self.active = True
        try:
            await asyncio.sleep(CALL)
            if self.fail_next:
                self.fail_next = False
                raise RuntimeError("worker failed")
            return answer
        finally:
            self.active = False

    async def summarize(self, text, backbone, require_mapped):
        return await self._call(json.dumps({"ok": True, "text": text}))

    async def setup(self, *args):
        return await self._call(json.dumps({"ok": False, "error": "stub"}))


class _Workers:
    def __init__(self):
        self.made = {}

    async def _get(self, name):
        return self.made.setdefault(name, _Worker())

    def __getitem__(self, name):
        return self._get(name)


def load_main():
    """Run main.py (minus its final `await start()`) against the stand-ins."""
    document, workers = _Document(), _Workers()
    stub = types.ModuleType("pyscript")
    stub.document, stub.workers = document, workers
    stub.when = lambda *args: (lambda func: func)
    source = MAIN.read_text()
    source = source[: source.rindex("await start()")]
    namespace = {"__name__": "main_under_test"}
    sys.modules["pyscript"] = stub
    try:
        exec(compile(source, str(MAIN), "exec"), namespace)
    finally:
        del sys.modules["pyscript"]
    namespace["LIVE_DEBOUNCE"] = DEBOUNCE
    return namespace, document, workers


class LiveCheckTests(unittest.IsolatedAsyncioTestCase):
    async def asyncSetUp(self):
        self.main, self.document, self.workers = load_main()
        await self.main["start"]()  # loads the workers and runs a first check
        self.live_worker = self.workers.made["w3"]  # the one the live check uses
        self.live_data = self.document.getElementById("live-data")
        self.live_data.textContent = ""

    def type(self, text):
        self.document.getElementById("input").value = text
        return asyncio.create_task(self.main["input_changed"](None))

    def shown(self):
        text = self.live_data.textContent
        return json.loads(text) if text else None

    async def until_with_worker(self):
        """Wait until a live check has a call out to its worker."""
        for _ in range(500):
            if self.live_worker.active:
                return
            await asyncio.sleep(0.01)
        self.fail("the live check never reached the worker")

    def assertNoOverlap(self):
        self.assertEqual(sum(w.overlaps for w in self.workers.made.values()), 0)

    async def test_edit_while_a_check_is_running(self):
        self.type("CD3,FITC")
        await self.until_with_worker()
        self.type("CD3,FITC\nCD4,APC")
        await asyncio.sleep(2.2)
        self.assertEqual(self.shown(), {"ok": True, "text": "CD3,FITC\nCD4,APC"})
        self.assertNoOverlap()

    async def test_burst_of_edits_ends_on_the_last_text(self):
        for i in range(20):
            self.type(f"line {i}")
            await asyncio.sleep(0.03)
        await asyncio.sleep(2.2)
        self.assertEqual(self.shown(), {"ok": True, "text": "line 19"})
        self.assertNoOverlap()

    async def test_compute_waits_for_a_running_check(self):
        self.type("CD3,FITC")
        await self.until_with_worker()
        self.main["run_btn"].disabled = False
        compute = asyncio.create_task(self.main["compute"](None))
        await asyncio.sleep(2.5)
        self.assertTrue(compute.done())
        self.assertNoOverlap()
        self.assertEqual(self.shown(), {"ok": True, "text": "CD3,FITC"})

    async def test_typing_during_compute_is_checked_afterwards(self):
        self.main["run_btn"].disabled = False
        asyncio.create_task(self.main["compute"](None))
        await asyncio.sleep(0.05)
        self.type("typed while computing")
        await asyncio.sleep(2.5)
        self.assertNoOverlap()
        self.assertEqual(self.shown(), {"ok": True, "text": "typed while computing"})

    async def test_worker_error_is_shown_and_next_edit_works(self):
        self.live_worker.fail_next = True
        self.type("A")
        await asyncio.sleep(1.2)
        self.assertFalse(self.shown()["ok"])
        self.assertIn("worker failed", self.shown()["error"])
        self.type("B")
        await asyncio.sleep(1.2)
        self.assertEqual(self.shown(), {"ok": True, "text": "B"})

    async def test_a_call_that_never_returns_is_given_up_on(self):
        self.main["LIVE_GIVE_UP"] = 1
        self.live_worker.hang_next = True
        self.type("A")
        await asyncio.sleep(1.7)  # the call has now been out longer than LIVE_GIVE_UP
        self.type("B")
        await asyncio.sleep(1.5)
        self.assertEqual(self.shown(), {"ok": True, "text": "B"})


if __name__ == "__main__":
    unittest.main()
