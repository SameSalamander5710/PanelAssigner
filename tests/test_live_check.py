"""The live check of the inputs (`refresh_live` in main.py) must never stay stale.

These tests load the real main.py with a stand-in for the `pyscript` module and
workers that answer one call at a time and never answer a call made while
another is out.

MicroPython in the browser runs a coroutine that JavaScript starts (an `async`
event handler, a top-level `await`) as a "top-level" coroutine, and all of those
share one wake-up slot: two suspended at once clobber each other (a call that is
never answered, an AttributeError about 'TopLevelCoro', a CPU spin). CPython's
asyncio does not do that, so it cannot be reproduced here; instead the rule that
avoids it is enforced: handlers are plain functions that return at once and the
module has no top-level `await` (see `entry_point_problems`). Handlers are called
the way JavaScript calls them, and the work they start runs as tasks.

Compute must not be held up by a live check: it starts at once, on the workers
that have no call out, and a worker that is still answering a check joins in
afterwards.
"""

import ast
import asyncio
import inspect
import json
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

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


SEARCH = 0.15  # how long one search task takes in a full compute
TASKS = 24


class _Worker:
    def __init__(self):
        self.active = False
        self.outstanding = 0  # calls not returned yet, including ones never answered
        self.overlaps = 0  # calls made while another call had not returned
        self.fail_next = False
        self.hang_next = False
        self.full = False  # setup/search/finalize behave like a real compute
        self.suggest_time = 0.01  # how long a suggest takes
        self.calls = []  # (name, time the call was made)
        self.never = []  # what calls that are never answered wait on; the browser keeps those alive

    async def _call(self, name, answer, duration=CALL):
        self.calls.append((name, asyncio.get_running_loop().time()))
        if self.outstanding:
            self.overlaps += 1
            self.outstanding += 1
            self.never.append(asyncio.Event())
            await self.never[-1].wait()  # never answered
        self.outstanding += 1
        try:
            if self.hang_next:
                self.hang_next = False
                self.never.append(asyncio.Event())
                await self.never[-1].wait()
            self.active = True
            await asyncio.sleep(duration)
            if self.fail_next:
                self.fail_next = False
                raise RuntimeError("worker failed")
            return answer
        finally:
            self.active = False
            self.outstanding -= 1

    def count(self, name):
        return sum(1 for call, _ in self.calls if call == name)

    async def summarize(self, text, backbone, require_mapped):
        return await self._call("summarize", json.dumps({"ok": True, "text": text}))

    async def setup(self, *args):
        if self.full:
            return await self._call("setup", json.dumps({"ok": True, "n_tasks": TASKS, "bound": "inf"}), 0.01)
        return await self._call("setup", json.dumps({"ok": False, "error": "stub"}))

    async def search(self, index, bound):
        return await self._call("search", "1.5\n{}", SEARCH)

    async def suggest(self, panel, n, exclude=""):
        return await self._call("suggest", json.dumps({"ok": True, "panel": panel, "n": n, "exclude": exclude}), self.suggest_time)

    async def finalize(self, results):
        return await self._call("finalize", json.dumps({
            "summary": "done", "notes": [], "html": "", "substitutions_html": "",
            "matrices_json": "{}"}), 0.01)


class _Workers:
    def __init__(self):
        self.made = {}

    async def _get(self, name):
        return self.made.setdefault(name, _Worker())

    def __getitem__(self, name):
        return self._get(name)


def entry_point_problems(source):
    """What would make MicroPython in the browser run a coroutine as a shared top-level one.

    Anything JavaScript starts must not suspend: a top-level `await` in the module, and
    `async` functions registered with @when, both do.
    """
    tree = ast.parse(source)
    problems = []

    def module_level_awaits(node):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda, ast.ClassDef)):
                continue
            if isinstance(child, ast.Await):
                problems.append(f"top-level await at line {child.lineno}")
            module_level_awaits(child)

    module_level_awaits(tree)
    for node in ast.walk(tree):
        if isinstance(node, ast.AsyncFunctionDef):
            for decorator in node.decorator_list:
                if isinstance(decorator, ast.Call) and getattr(decorator.func, "id", "") == "when":
                    problems.append(f"async @when handler {node.name} at line {node.lineno}")
    return problems


def load_main():
    """Run main.py against the stand-ins (including its final `spawn(start())`)."""
    document, workers = _Document(), _Workers()
    stub = types.ModuleType("pyscript")
    stub.document, stub.workers = document, workers
    stub.when = lambda *args: (lambda func: func)
    source = MAIN.read_text()
    problems = entry_point_problems(source)
    if problems:
        raise AssertionError("main.py would break MicroPython in the browser: " + "; ".join(problems))
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
        self.tasks = []  # every task main.py starts; one nobody refers to can be collected while it waits
        real_create_task = asyncio.create_task

        def create_task(coro):
            task = real_create_task(coro)
            self.tasks.append(task)
            return task

        patcher = mock.patch.object(asyncio, "create_task", create_task)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.main, self.document, self.workers = load_main()
        await self.idle()  # start() loaded the workers and ran a first check
        self.live_worker = self.workers.made["w3"]  # the one the live check uses
        self.live_data = self.document.getElementById("live-data")
        self.live_data.textContent = ""

    async def idle(self):
        await self.tasks[0]  # start()
        while self.main["live"]["running"]:
            await asyncio.sleep(0.01)

    def js(self, handler, event=None):
        """Call an event handler the way JavaScript does; it must be done at once."""
        result = self.main[handler](event)
        if inspect.iscoroutine(result):
            result.close()
        self.assertIsNone(result, handler + " must return at once, not a coroutine (see entry_point_problems)")

    def type(self, text):
        self.document.getElementById("input").value = text
        self.js("input_changed")

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

    def start_full_compute(self):
        """Click Compute with workers that really search; returns the task that computes."""
        for worker in self.workers.made.values():
            worker.full = True
        self.main["run_btn"].disabled = False
        before = len(self.tasks)
        self.js("compute")
        self.assertEqual(len(self.tasks), before + 1)  # the click starts one task
        return self.tasks[-1]

    async def test_compute_does_not_wait_for_a_running_check(self):
        self.type("CD3,FITC")
        await self.until_with_worker()  # the check is with w3 for the next CALL seconds
        began = asyncio.get_running_loop().time()
        compute = self.start_full_compute()
        await asyncio.sleep(0.1)
        started = [t for w in self.workers.made.values() for name, t in w.calls if name == "setup"]
        self.assertTrue(started, "Compute has not reached a worker")
        self.assertLess(min(started) - began, 0.05)
        self.assertEqual(self.live_worker.count("setup"), 0)  # w3 is not sent anything yet
        await compute
        self.assertNoOverlap()

    async def test_the_busy_worker_joins_once_its_check_has_returned(self):
        self.type("CD3,FITC")
        await self.until_with_worker()
        await self.start_full_compute()
        self.assertEqual(self.live_worker.count("setup"), 1)
        self.assertGreater(self.live_worker.count("search"), 0)
        self.assertNoOverlap()

    async def test_check_answer_is_dropped_during_compute_and_summary_refreshed_after(self):
        self.type("CD3,FITC")
        await self.until_with_worker()
        compute = self.start_full_compute()
        await asyncio.sleep(CALL + 0.2)  # the check has returned, Compute is still going
        self.assertFalse(compute.done())
        self.assertIsNone(self.shown())  # nothing is shown while Compute runs
        await compute
        await asyncio.sleep(CALL + 0.4)
        self.assertEqual(self.shown(), {"ok": True, "text": "CD3,FITC"})
        self.assertNoOverlap()

    async def test_click_while_waiting_for_typing_to_pause(self):
        self.type("CD3,FITC")  # the check is still in its pause
        began = asyncio.get_running_loop().time()
        compute = self.start_full_compute()
        await asyncio.sleep(0.05)
        for name in ("w0", "w1", "w2", "w3"):  # no check was started, so all are free
            self.assertEqual(self.workers.made[name].count("setup"), 1, name)
        self.assertLess(self.workers.made["w3"].calls[-1][1] - began, 0.1)
        await compute
        self.assertEqual(self.workers.made["w3"].count("summarize") >= 0, True)
        await asyncio.sleep(CALL + 0.4)
        self.assertEqual(self.shown(), {"ok": True, "text": "CD3,FITC"})
        self.assertNoOverlap()

    async def test_compute_does_not_wait_for_a_check_that_never_returns(self):
        self.live_worker.hang_next = True
        self.type("CD3,FITC")
        await asyncio.sleep(DEBOUNCE + 0.1)  # the check is now out with w3 and will never return
        began = asyncio.get_running_loop().time()
        await self.start_full_compute()
        self.assertLess(asyncio.get_running_loop().time() - began, 2.5)  # not LIVE_GIVE_UP
        self.assertEqual(self.live_worker.count("setup"), 0)
        self.assertGreater(self.workers.made["w1"].count("search"), 0)
        await asyncio.sleep(CALL + 0.4)  # the summary comes from another worker
        self.assertEqual(self.shown(), {"ok": True, "text": "CD3,FITC"})
        self.assertNoOverlap()

    async def test_typing_during_compute_is_checked_afterwards(self):
        self.main["run_btn"].disabled = False
        self.js("compute")
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

    async def test_suggest_returns_at_once_and_answers_from_its_task(self):
        self.document.getElementById("suggest-panel").textContent = "2"
        self.document.getElementById("suggest-n").value = "1"
        self.js("suggest")
        self.assertTrue(self.document.getElementById("suggest-run").disabled)
        await asyncio.sleep(0.2)
        self.assertEqual(json.loads(self.document.getElementById("suggest-data").textContent)["panel"], 2)
        self.assertFalse(self.document.getElementById("suggest-run").disabled)

    async def test_suggest_passes_the_leave_out_box_to_the_worker(self):
        self.document.getElementById("suggest-panel").textContent = "1"
        self.document.getElementById("suggest-n").value = "2"
        self.document.getElementById("suggest-exclude").value = "BV*\nQdot 605"
        self.js("suggest")
        await asyncio.sleep(0.2)
        answer = json.loads(self.document.getElementById("suggest-data").textContent)
        self.assertEqual(answer["exclude"], "BV*\nQdot 605")

    def start_suggest(self, duration):
        """Click Suggest with worker 0 taking `duration` seconds to answer."""
        self.workers.made["w0"].suggest_time = duration
        self.document.getElementById("suggest-panel").textContent = "1"
        self.document.getElementById("suggest-n").value = "3"
        self.js("suggest")

    async def test_compute_is_not_started_while_a_suggest_is_out(self):
        # Worker 0 holds the panels. A Compute that left it out of the search (it was busy)
        # and then sent it finalize worked from the previous run's data, on top of the
        # Suggest that was still running there; the call was never answered.
        w0 = self.workers.made["w0"]
        self.start_suggest(0.6)
        await asyncio.sleep(0.05)
        self.assertTrue(self.main["run_btn"].disabled)  # Compute cannot be clicked meanwhile
        for worker in self.workers.made.values():
            worker.full = True
        before = len(self.tasks)
        self.js("compute")  # a click that gets through anyway does nothing
        self.assertEqual(len(self.tasks), before)
        await asyncio.sleep(0.8)
        self.assertFalse(self.main["run_btn"].disabled)
        self.assertFalse(self.document.getElementById("suggest-run").disabled)
        self.assertEqual([name for name, _ in w0.calls], ["suggest"])
        self.assertNoOverlap()

    async def test_compute_after_a_suggest_sets_up_worker_zero_before_it_finalizes(self):
        w0 = self.workers.made["w0"]
        self.start_suggest(0.3)
        await asyncio.sleep(0.5)
        await self.start_full_compute()
        names = [name for name, _ in w0.calls]
        self.assertEqual(names[0], "suggest")
        self.assertIn("setup", names)
        self.assertLess(names.index("setup"), names.index("finalize"))
        self.assertEqual(self.document.getElementById("status").innerText[:4], "done")
        self.assertNoOverlap()

    async def test_suggest_is_not_started_while_compute_runs(self):
        w0 = self.workers.made["w0"]
        compute = self.start_full_compute()
        await asyncio.sleep(0.1)
        self.assertTrue(self.document.getElementById("suggest-run").disabled)
        self.js("suggest")  # a click that gets through anyway does nothing
        self.assertEqual(w0.count("suggest"), 0)
        self.assertIn("Compute", self.document.getElementById("suggest-status").innerText)
        await compute
        self.assertFalse(self.document.getElementById("suggest-run").disabled)
        self.start_suggest(0.01)
        await asyncio.sleep(0.2)
        self.assertEqual(w0.count("suggest"), 1)
        self.assertNoOverlap()

    async def test_a_failed_suggest_gives_both_buttons_back(self):
        self.workers.made["w0"].fail_next = True
        self.start_suggest(0.01)
        await asyncio.sleep(CALL + 0.3)
        self.assertIn("worker failed", self.document.getElementById("suggest-status").innerText)
        self.assertFalse(self.document.getElementById("suggest-run").disabled)
        self.assertFalse(self.main["run_btn"].disabled)

    async def test_suggest_returns_the_live_check_that_it_paused(self):
        self.start_suggest(0.05)
        self.type("typed while suggesting")  # run_btn is disabled, so no check is made
        await asyncio.sleep(0.4)
        self.assertIsNone(self.shown())
        await asyncio.sleep(CALL + 0.4)  # ... the check runs once Suggest is done
        self.assertEqual(self.shown(), {"ok": True, "text": "typed while suggesting"})
        self.assertNoOverlap()

    async def test_every_handler_returns_at_once(self):
        for handler in ("input_changed", "backbone_changed", "require_mapped_changed"):
            self.js(handler)
        await asyncio.sleep(0.2)


class EntryPointRuleTests(unittest.TestCase):
    def test_main_py_follows_the_rule(self):
        self.assertEqual(entry_point_problems(MAIN.read_text()), [])

    def test_the_check_catches_what_broke_the_page(self):
        self.assertEqual(entry_point_problems("async def helper():\n    await x()\n"), [])
        self.assertEqual(len(entry_point_problems("await start()\n")), 1)
        self.assertEqual(len(entry_point_problems(
            "@when('input', '#a')\nasync def changed(event):\n    await refresh()\n")), 1)
        self.assertEqual(len(entry_point_problems(
            "@when('input', '#a')\ndef changed(event):\n    spawn(refresh())\n")), 0)


if __name__ == "__main__":
    unittest.main()
