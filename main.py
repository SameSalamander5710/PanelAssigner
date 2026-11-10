import asyncio
import json
import time
from pyscript import document, workers, when

NAMES = ["w0", "w1", "w2", "w3"]  # must match the <script ... name="..."> tags
pool = []
status = document.getElementById("status")
run_btn = document.getElementById("run")


async def start():
    for name in NAMES:
        pool.append(await workers[name])  # resolves once that worker has loaded
    run_btn.disabled = False
    status.innerText = "Ready (" + str(len(pool)) + " workers)."


async def drive(worker, text, top_n, state):  # top_n is an int here
    """One loop per worker: set up, then keep pulling the next unclaimed subtree."""
    info = json.loads(await worker.setup(text, str(top_n)))
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
    out.value = ""
    status.innerText = "Computing on " + str(len(pool)) + " workers..."
    t0 = time.time()
    state = {"next": 0, "results": [], "known": [], "bound": None, "error": None}
    text = document.getElementById("input").value
    top_n = int(document.getElementById("topn").value)
    await asyncio.gather(*[drive(w, text, top_n, state) for w in pool])

    if state["error"]:
        status.innerText = "Error: " + state["error"]
    else:
        res = json.loads(await pool[0].finalize("[" + ",".join(state["results"]) + "]"))
        status.innerText = res["summary"] + " Time: " + str(round(time.time() - t0, 1)) + " s."
        out.value = res["tsv"]
    run_btn.disabled = False


await start()
