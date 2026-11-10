"""One of several identical PyScript (Pyodide) workers. Each runs in its own
Web Worker thread; main.py hands out search subtrees to whichever is free."""
import csv
import html as html_lib
import io
import json

import Panel_optimized as P

SIM_DF, SIM_NAMES = P.load_similarity_matrix("./constants/Similarity_Index.csv")
DET_MAP = P.load_detector_map("./constants/Detector_Fluorophores.csv")

TARGET_TASKS = 32  # more, smaller subtrees = better load balancing and pruning
STATE = {}


def parse_pasted(text):
    """Parse pasted Antibody/Fluorophore pairs (tab, comma or semicolon separated)."""
    text = text.strip()
    delim = "\t" if "\t" in text else ("," if "," in text else ";")
    antibodies = {}
    for row in csv.reader(io.StringIO(text), delimiter=delim):
        if len(row) < 2:
            continue
        ab, fl = row[0].strip(), row[1].strip()
        if not ab or not fl:
            continue
        if ab.lower() == "antibody" and fl.lower() == "fluorophore":
            continue
        options = antibodies.setdefault(ab, [])
        if fl not in options:
            options.append(fl)
    if not antibodies:
        raise ValueError("No Antibody/Fluorophore pairs found. Paste two columns.")
    return antibodies


def setup(text, top_n=3):
    """Every worker runs this with the same input and builds the same task list."""
    try:
        top_n = int(top_n)
        antibodies = parse_pasted(text)
        data = P.prepare_search(antibodies, SIM_DF, SIM_NAMES, DET_MAP, top_n)
        if data is None:
            raise ValueError("No assignable antibody/fluorophore pairs.")
        tasks = P.make_search_tasks(data, TARGET_TASKS)
        STATE.update(data=data, tasks=tasks, top_n=top_n, total=len(antibodies),
                     max_count=data["max_count"])
        return json.dumps({"ok": True, "n_tasks": len(tasks),
                           "bound": repr(P.seed_bound(data))})
    except Exception as exc:
        return json.dumps({"ok": False, "error": str(exc)})


def search(i, bound="inf"):
    """Search one subtree under the current global bound.

    Returns "c1,c2,...\n<json>": the first line lists this subtree's best
    complexities (so the page can tighten the bound for the other workers).
    """
    b = float(bound)
    if b != float("inf"):
        b *= 1 + 1e-6  # the bound travelled as text; never let rounding prune a tie
    best, evaluated, pruned = P.search_subtree(STATE["data"], STATE["tasks"][i], b)
    out = []
    for complexity, fluor_set, panel in best:
        panel = dict(panel)
        panel["unmapped_fluors"] = sorted(panel["unmapped_fluors"])
        out.append([complexity, list(fluor_set), panel])
    head = ",".join(repr(c) for c, _, _ in best)
    return head + "\n" + json.dumps([out, evaluated, pruned])


def finalize(results_json):
    """Merge all subtree results and format them (run on one worker only)."""
    results = json.loads(results_json)
    evaluated = sum(r[1] for r in results)
    pruned = sum(r[2] for r in results)
    best = [(c, tuple(fs), panel) for local, _, _ in results for c, fs, panel in local]
    best.sort(key=lambda x: (x[0], x[1]))
    best = best[: STATE["top_n"]]

    def esc(value):
        return html_lib.escape(str(value), quote=False)

    parts = [
        '<table class="results-table"><thead><tr>',
        '<th class="num">#</th><th>Antibody</th>',
    ]
    for i in range(1, len(best) + 1):
        parts.append(f"<th>Panel {i}</th>")
    parts.append("</tr></thead><tbody>")

    for row_no, ab in enumerate(STATE["data"]["original_order"], 1):
        parts.append(f'<tr><td class="num">{row_no}</td><td>{esc(ab)}</td>')
        for _, _, panel in best:
            if ab in panel["assignments"]:
                fl = panel["assignments"][ab]
                if fl in panel["unmapped_fluors"]:
                    cls, title = "unmapped", ' title="unmapped"'
                else:
                    cls, title = "assigned", ""
                    if fl in panel["substitutions"]:
                        sub = html_lib.escape(panel["substitutions"][fl])
                        title = f' title="substituted: {sub}"'
                parts.append(f'<td class="{cls}"{title}>{esc(fl)}</td>')
            else:
                parts.append('<td class="unassigned" title="unassigned"></td>')
        parts.append("</tr>")

    parts.append('</tbody><tfoot><tr class="ci-row"><th colspan="2">Matrix_CI</th>')
    for c, _, _ in best:
        parts.append(f'<td class="num">{c:.4f}</td>')
    parts.append("</tr></tfoot></table>")

    summary = (
        f"Maximum panel size: {STATE['max_count']}/{STATE['total']}. "
        f"Fluorophore sets scored: {evaluated:,} (branches pruned: {pruned:,})."
    )
    return json.dumps({"summary": summary, "html": "".join(parts)})


__export__ = ["setup", "search", "finalize"]
