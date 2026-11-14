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


def set_mapped_info(notes, has_unmapped=False):
    """Put the per-run notes (unmapped, dropped, merged) in the hover text."""
    base = mapped_info.getAttribute("data-info") or ""
    mapped_tip.textContent = base + ("\n\n" + "\n".join(notes) if notes else "")
    if has_unmapped:
        mapped_info.classList.add("has-unmapped")
    else:
        mapped_info.classList.remove("has-unmapped")


async def refresh_live():
    """Check the inputs on a worker once typing pauses; index.html shows the answer.

    At most one check loop runs, and it has at most one call out to the worker. Editing
    while a call is running must not start a second call: the worker answers one call
    at a time, and an overlapping call was never answered, which left the summary
    greyed out for good. Instead the loop notices the newer text when the call returns,
    and checks again (after the usual pause).
    """
    live["wanted"] += 1
    if live["running"]:
        if time.time() - live["since"] < LIVE_GIVE_UP:
            return  # the running loop will pick up this edit
        # It has made no progress for a long time: abandon it and start afresh.
    live["gen"] += 1
    mine = live["gen"]
    live["running"] = True
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
            raw = await pool[-1].summarize(  # the last worker: w0 holds the panels for Suggest
                document.getElementById("input").value,
                document.getElementById("backbone").value,
                "1" if document.getElementById("requiremapped").checked else "",
            )
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


async def live_check_idle():
    """Wait for a running live check to finish its call to the worker, so that
    nothing else is sent to that worker at the same time."""
    give_up = time.time() + LIVE_GIVE_UP
    while live["running"] and time.time() < give_up:
        await asyncio.sleep(0.05)


async def start():
    for name in NAMES:
        pool.append(await workers[name])  # resolves once that worker has loaded
    run_btn.disabled = False
    status.innerText = "Ready (" + str(len(pool)) + " workers)."
    await refresh_live()  # the browser may have restored text from an earlier visit


async def drive(worker, text, backbone, top_n, max_fit, require_mapped, metric, state):  # top_n is an int here
    """One loop per worker: set up, then keep pulling the next unclaimed subtree."""
    info = json.loads(await worker.setup(text, str(top_n), backbone, max_fit, require_mapped, metric))
    if not info["ok"]:
        state["error"] = info["error"]
        return
    if state["bound"] is None:
        state["bound"] = float(info["bound"])
    n = info["n_tasks"]
    while state["next"] < n and not state["error"]:
        i = state["next"]
        state["next"] = i + 1
        raw = await worker.search(i, str(state["bound"]))
        head, _, body = raw.partition("\n")
        state["results"].append(body)
        # tighten the shared bound so the other workers can prune harder
        for s in head.split(","):
            if s:
                state["known"].append(float(s))
        state["known"].sort()
        if len(state["known"]) >= top_n:
            kth = state["known"][top_n - 1]
            if kth < state["bound"]:
                state["bound"] = kth


@when("input", "#input")
async def input_changed(event):
    await refresh_live()


@when("input", "#backbone")
async def backbone_changed(event):
    await refresh_live()


@when("change", "#requiremapped")
async def require_mapped_changed(event):
    await refresh_live()


@when("click", "#run")
async def compute(event):
    run_btn.disabled = True  # also ends a live check that is waiting for typing to pause
    await live_check_idle()  # one that is already with a worker finishes first
    out = document.getElementById("output")
    out.innerHTML = ""
    subs.innerHTML = ""
    matrix_data.textContent = ""
    suggest_data.textContent = ""
    suggest_status.innerText = ""
    set_mapped_info([])
    status.innerText = "Computing on " + str(len(pool)) + " workers..."
    t0 = time.time()
    state = {"next": 0, "results": [], "known": [], "bound": None, "error": None}
    text = document.getElementById("input").value
    backbone = document.getElementById("backbone").value
    top_n = int(document.getElementById("topn").value)
    max_fit = document.getElementById("maxfit").value  # "" = no limit
    require_mapped = "1" if document.getElementById("requiremapped").checked else ""
    metric = document.getElementById("metric").value  # key of Panel_optimized.METRICS
    await asyncio.gather(*[drive(w, text, backbone, top_n, max_fit, require_mapped, metric, state) for w in pool])

    if state["error"]:
        status.innerText = "Error: " + state["error"]
    else:
        res = json.loads(await pool[0].finalize("[" + ",".join(state["results"]) + "]"))
        status.innerText = res["summary"] + " Time: " + str(round(time.time() - t0, 1)) + " s."
        set_mapped_info(res.get("notes") or [], res.get("has_unmapped", False))
        out.innerHTML = res["html"]
        subs.innerHTML = res["substitutions_html"]
        matrix_data.textContent = res["matrices_json"]
    run_btn.disabled = False
    await refresh_live()  # the text may have been edited while computing


@when("click", "#suggest-run")
async def suggest(event):
    suggest_btn.disabled = True
    suggest_data.textContent = ""
    try:
        panel_no = int(document.getElementById("suggest-panel").textContent)
        n = int(document.getElementById("suggest-n").value)
    except ValueError:
        suggest_status.innerText = "Enter how many fluorophores to add (1-5)."
        suggest_btn.disabled = False
        return
    suggest_status.innerText = "Searching..."
    raw = await pool[0].suggest(panel_no, n)  # same worker that holds the panels
    res = json.loads(raw)
    if res["ok"]:
        suggest_status.innerText = ""
        suggest_data.textContent = raw
    else:
        suggest_status.innerText = "Error: " + res["error"]
    suggest_btn.disabled = False


await start()
