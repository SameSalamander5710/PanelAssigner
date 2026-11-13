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


async def start():
    for name in NAMES:
        pool.append(await workers[name])  # resolves once that worker has loaded
    run_btn.disabled = False
    status.innerText = "Ready (" + str(len(pool)) + " workers)."


async def drive(worker, text, backbone, top_n, max_fit, require_mapped, state):  # top_n is an int here
    """One loop per worker: set up, then keep pulling the next unclaimed subtree."""
    info = json.loads(await worker.setup(text, str(top_n), backbone, max_fit, require_mapped))
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


@when("click", "#run")
async def compute(event):
    run_btn.disabled = True
    out = document.getElementById("output")
    out.innerHTML = ""
    subs.innerHTML = ""
    matrix_data.textContent = ""
    suggest_data.textContent = ""
    suggest_status.innerText = ""
    status.innerText = "Computing on " + str(len(pool)) + " workers..."
    t0 = time.time()
    state = {"next": 0, "results": [], "known": [], "bound": None, "error": None}
    text = document.getElementById("input").value
    backbone = document.getElementById("backbone").value
    top_n = int(document.getElementById("topn").value)
    max_fit = document.getElementById("maxfit").value  # "" = no limit
    require_mapped = "1" if document.getElementById("requiremapped").checked else ""
    await asyncio.gather(*[drive(w, text, backbone, top_n, max_fit, require_mapped, state) for w in pool])

    if state["error"]:
        status.innerText = "Error: " + state["error"]
    else:
        res = json.loads(await pool[0].finalize("[" + ",".join(state["results"]) + "]"))
        status.innerText = res["summary"] + " Time: " + str(round(time.time() - t0, 1)) + " s."
        out.innerHTML = res["html"]
        subs.innerHTML = res["substitutions_html"]
        matrix_data.textContent = res["matrices_json"]
    run_btn.disabled = False


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
