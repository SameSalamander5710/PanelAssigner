import asyncio
import json
import time
from pyscript import document, workers, when

NAMES = ["w0", "w1", "w2", "w3"]  # must match the <script ... name="..."> tags
pool = []
status = document.getElementById("status")
run_btn = document.getElementById("run")
subs = document.getElementById("substitutions")
matrix_data = document.getElementById("matrix-data")  # index.html draws the plot from this
suggest_btn = document.getElementById("suggest-run")
suggest_status = document.getElementById("suggest-status")
suggest_data = document.getElementById("suggest-data")  # ... and the suggestions from this
mapped_info = document.getElementById("requiremapped-info")  # hover notes for the matrix option
mapped_tip = document.getElementById("requiremapped-tip")
live_data = document.getElementById("live-data")  # ... and the live summary of the inputs from this
LIVE_DEBOUNCE = 0.3  # seconds without typing before the inputs are checked
LIVE_GIVE_UP = 15  # seconds; a check that takes longer than this is abandoned
# wanted: counts edits. running: a check loop is active. gen: which loop is current.
# since: when the loop last made progress.
live = {"wanted": 0, "running": False, "gen": 0, "since": 0}
# Calls a worker has out. Nothing is sent to a worker that has a call out (overlapping
# calls were once never answered; see spawn() for the cause that was found later, and
# this is kept as a precaution). So Compute leaves such a worker to join in later
# instead of waiting for it.
busy = {}
# Worker 0 holds the panels that Suggest reads, and Compute rebuilds them there (setup,
# then finalize). So only one of the two runs at a time: the button of the other is
# disabled meanwhile. Without this, a Compute started during a slow Suggest left worker 0
# out of the search (it was busy), then sent it finalize with the previous run's data, on
# top of the Suggest still running there.
active = {"compute": False, "suggest": False}


def spawn(coro):
    """Run `coro` as a task of its own.

    MicroPython in the browser runs a coroutine that JavaScript starts (an `async` event
    handler, or a top-level `await`) as a "top-level" coroutine, and all of those share
    ONE wake-up slot. If one is suspended (in a sleep, a worker call, gather) and another
    starts and suspends, they overwrite each other: one is never resumed (a summary that
    never arrives, a Compute that hangs), the scheduler throws into the shared slot
    ("TopLevelCoro has no attribute 'throw'"), or the CPU spins.

    So nothing that JavaScript starts may suspend. Event handlers and the module itself
    do only synchronous work and hand the rest to spawn(); a Task has its own identity.
    tests/test_live_check.py checks that no handler is `async` and that the module has
    no top-level `await`.
    """
    return asyncio.create_task(coro)


async def call(i, function, *args):
    """Call `function` (a method of worker i) and count it as out until it returns."""
    busy[i] = busy.get(i, 0) + 1
    try:
        return await function(*args)
    finally:
        busy[i] -= 1


def live_worker():
    """A free worker for the live check: the last one if it is free. Never the first,
    which holds the panels for Suggest."""
    for i in range(len(pool) - 1, 0, -1):
        if not busy.get(i):
            return i
    return None


def set_mapped_info(notes, has_unmapped=False):
    """Put the per-run notes (unmapped, dropped, merged) in the hover text."""
    base = mapped_info.getAttribute("data-info") or ""
    mapped_tip.textContent = base + ("\n\n" + "\n".join(notes) if notes else "")
    if has_unmapped:
        mapped_info.classList.add("has-unmapped")
    else:
        mapped_info.classList.remove("has-unmapped")


def refresh_live(fresh=False):
    """Check the inputs on a worker once typing pauses; index.html shows the answer.

    Returns at once: the check runs as a task (see spawn()). At most one check loop
    runs, and it has at most one call out to a worker. Editing while a call is running
    does not start a second call; the loop notices the newer text when the call returns,
    and checks again (after the usual pause).

    `fresh` replaces a running loop (run_compute() does this when it is done, so that a
    check that was left with a worker meanwhile cannot hold up the summary).
    """
    live["wanted"] += 1
    if live["running"] and not fresh:
        if time.time() - live["since"] < LIVE_GIVE_UP:
            return  # the running loop will pick up this edit
        # It has made no progress for a long time: abandon it and start afresh.
    live["gen"] += 1
    live["running"] = True
    spawn(live_loop(live["gen"]))


async def live_loop(mine):
    try:
        while live["gen"] == mine:
            seen = live["wanted"]
            live["since"] = time.time()
            await asyncio.sleep(LIVE_DEBOUNCE)
            if live["wanted"] != seen:
                continue  # typing again; wait for the next pause
            if not pool or run_btn.disabled:  # not ready, or a compute is running
                break  # start() and compute() ask again when they are done
            live["since"] = time.time()
            i = live_worker()
            if i is None:
                raise RuntimeError("no free worker")
            raw = await call(
                i, pool[i].summarize,
                document.getElementById("input").value,
                document.getElementById("backbone").value,
                "1" if document.getElementById("requiremapped").checked else "",
            )
            if run_btn.disabled:  # Compute was clicked meanwhile; it refreshes the summary when done
                break
            if live["gen"] != mine:
                break  # abandoned while waiting
            if live["wanted"] == seen:  # nothing was typed meanwhile
                live_data.textContent = raw
                break
            # Edited during the call: go round again, with the new text.
    except Exception as exc:
        if live["gen"] == mine:  # shown under the inputs instead of staying greyed out
            live_data.textContent = json.dumps(
                {"ok": False, "error": "Could not check the input: " + str(exc)})
    finally:
        if live["gen"] == mine:
            live["running"] = False


async def start():
    for name in NAMES:
        pool.append(await workers[name])  # resolves once that worker has loaded
    run_btn.disabled = False
    status.innerText = "Ready (" + str(len(pool)) + " workers)."
    refresh_live()  # the browser may have restored text from an earlier visit


async def drive(i, text, backbone, top_n, max_fit, require_mapped, metric, state):  # top_n is an int here
    """One loop per worker: set up, then keep pulling the next unclaimed subtree."""
    worker = pool[i]
    info = json.loads(await call(i, worker.setup, text, str(top_n), backbone, max_fit, require_mapped, metric))
    if not info["ok"]:
        state["error"] = info["error"]
        return
    if state["bound"] is None:
        state["bound"] = float(info["bound"])
    n = state["n"] = info["n_tasks"]
    while state["next"] < n and not state["error"]:
        k = state["next"]
        state["next"] = k + 1
        raw = await call(i, worker.search, k, str(state["bound"]))
        head, _, body = raw.partition("\n")
        state["results"].append(body)
        # tighten the shared bound so the other workers can prune harder
        for v in head.split(","):
            if v:
                state["known"].append(float(v))
        state["known"].sort()
        if len(state["known"]) >= top_n:
            kth = state["known"][top_n - 1]
            if kth < state["bound"]:
                state["bound"] = kth


async def drive_now(i, *args):
    """drive() on a worker that was free when Compute was clicked."""
    state = args[-1]
    try:
        await drive(i, *args)
    finally:
        state["active"] -= 1


async def drive_later(i, *args):
    """A worker that still had a live check out when Compute was clicked. It joins in
    once it has answered, which is normally a moment. If the others finish first, or it
    never answers, Compute goes on without it."""
    state = args[-1]
    while busy.get(i) and state["active"] > 0:
        await asyncio.sleep(0.05)
    if busy.get(i) or state["error"] or (state["n"] is not None and state["next"] >= state["n"]):
        return
    await drive(i, *args)


# Event handlers are plain functions that return at once (see spawn()).
@when("input", "#input")
def input_changed(event):
    refresh_live()


@when("input", "#backbone")
def backbone_changed(event):
    refresh_live()


@when("change", "#requiremapped")
def require_mapped_changed(event):
    refresh_live()


@when("click", "#run")
def compute(event):
    # Compute starts at once. A live check that is waiting for typing to pause ends now
    # (run_btn is disabled when it wakes), and the answer of one that is already with a
    # worker is dropped; the summary is refreshed when Compute is done.
    if active["suggest"]:  # the button is disabled meanwhile; this is only a guard
        return
    active["compute"] = True
    run_btn.disabled = True
    suggest_btn.disabled = True  # Suggest reads what this Compute is about to replace
    spawn(run_compute())


async def run_compute():
    out = document.getElementById("output")
    out.innerHTML = ""
    subs.innerHTML = ""
    matrix_data.textContent = ""
    suggest_data.textContent = ""
    suggest_status.innerText = ""
    set_mapped_info([])
    status.innerText = "Computing on " + str(len(pool)) + " workers..."
    t0 = time.time()
    state = {"next": 0, "n": None, "active": 0, "results": [], "known": [], "bound": None, "error": None}
    text = document.getElementById("input").value
    backbone = document.getElementById("backbone").value
    top_n = int(document.getElementById("topn").value)
    max_fit = document.getElementById("maxfit").value  # "" = no limit
    require_mapped = "1" if document.getElementById("requiremapped").checked else ""
    metric = document.getElementById("metric").value  # key of Panel_optimized.METRICS
    free = [i for i in range(len(pool)) if not busy.get(i)]
    later = [i for i in range(len(pool)) if busy.get(i)]
    if not free:  # nothing is free (not expected): use them all
        free, later = later, []
    state["active"] = len(free)
    args = (text, backbone, top_n, max_fit, require_mapped, metric, state)
    await asyncio.gather(*([drive_now(i, *args) for i in free] + [drive_later(i, *args) for i in later]))

    if state["error"]:
        status.innerText = "Error: " + state["error"]
    else:
        res = json.loads(await call(0, pool[0].finalize, "[" + ",".join(state["results"]) + "]"))
        status.innerText = res["summary"] + " Time: " + str(round(time.time() - t0, 1)) + " s."
        set_mapped_info(res.get("notes") or [], res.get("has_unmapped", False))
        out.innerHTML = res["html"]
        subs.innerHTML = res["substitutions_html"]
        matrix_data.textContent = res["matrices_json"]
    active["compute"] = False
    run_btn.disabled = False
    suggest_btn.disabled = False
    refresh_live(True)  # the text may have been edited while computing; and this is the refresh after Compute


@when("click", "#suggest-run")
def suggest(event):
    if not pool or active["compute"] or active["suggest"]:
        if active["compute"]:
            suggest_status.innerText = "Wait for Compute to finish."
        return  # the button is disabled meanwhile; this is only a guard
    suggest_btn.disabled = True
    suggest_data.textContent = ""
    try:
        panel_no = int(document.getElementById("suggest-panel").textContent)
        n = int(document.getElementById("suggest-n").value)
    except ValueError:
        suggest_status.innerText = "Enter how many fluorophores to add (1-5)."
        suggest_btn.disabled = False
        return
    exclude = document.getElementById("suggest-exclude").value  # fluorophores to leave out
    active["suggest"] = True
    run_btn.disabled = True  # Compute would take worker 0 from under this Suggest
    suggest_status.innerText = "Searching... (Compute is available again when this is done)"
    spawn(run_suggest(panel_no, n, exclude))


async def run_suggest(panel_no, n, exclude):
    try:
        raw = await call(0, pool[0].suggest, panel_no, n, exclude)  # same worker that holds the panels
        res = json.loads(raw)
        if res["ok"]:
            suggest_status.innerText = ""
            suggest_data.textContent = raw
        else:
            suggest_status.innerText = "Error: " + res["error"]
    except Exception as exc:  # shown here instead of leaving both buttons disabled
        suggest_status.innerText = "Error: " + str(exc)
    finally:
        active["suggest"] = False
        suggest_btn.disabled = False
        run_btn.disabled = False
        refresh_live(True)  # a check that was dropped while Suggest held Compute's button


spawn(start())
