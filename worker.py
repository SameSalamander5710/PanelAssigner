"""One of several identical PyScript (Pyodide) workers. Each runs in its own
Web Worker thread; main.py hands out search subtrees to whichever is free."""
import csv
import html as html_lib
import io
import json
import math

import Panel_optimized as P

SIM_DF, SIM_NAMES = P.load_similarity_matrix("./constants/Similarity_Index.csv")
DET_MAP = P.load_detector_map("./constants/Detector_Fluorophores.csv")
try:
    ALIASES, ALIAS_ERROR = P.load_fluorophore_aliases("./constants/Fluorophore_Aliases.csv"), None
except Exception as exc:  # a typo in the CSV is reported in the UI, not lost in the worker
    ALIASES, ALIAS_ERROR = [], str(exc)

TARGET_TASKS = 32  # more, smaller subtrees = better load balancing and pruning
STATE = {}


def parse_pasted(text, required=True):
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
    if not antibodies and required:
        raise ValueError("No Antibody/Fluorophore pairs found. Paste two columns.")
    return antibodies


def parse_backbone(text):
    """Parse optional fixed backbone pairs; the first row per antibody wins."""
    text = (text or "").strip()
    if not text:
        return {}
    pairs = parse_pasted(text, required=False)
    if not pairs:
        raise ValueError("No Antibody/Fluorophore pairs found in the fixed backbone.")
    return {ab: fls[0] for ab, fls in pairs.items() if fls}


def setup(text, top_n=3, backbone_text=""):
    """Every worker runs this with the same input and builds the same task list."""
    try:
        if ALIAS_ERROR:
            raise ValueError(f"Could not load fluorophore aliases: {ALIAS_ERROR}")
        top_n = int(top_n)
        backbone = parse_backbone(backbone_text)
        antibodies = parse_pasted(text, required=not backbone)
        data = P.prepare_search(
            antibodies, SIM_DF, SIM_NAMES, DET_MAP, top_n,
            backbone=backbone, aliases=ALIASES,
        )
        if data is None:
            raise ValueError("No assignable antibody/fluorophore pairs.")
        tasks = P.make_search_tasks(data, TARGET_TASKS)
        STATE.update(data=data, tasks=tasks, top_n=top_n,
                     total=len(data["original_order"]),
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

    # One similarity matrix per panel for the plot below the results. Names go
    # out as JSON; the page inserts them as text, never as HTML.
    matrices = []
    for complexity, _, panel in best:
        labels = list(panel["assignments"].values())  # same order as the table
        backbone_fluors = set(panel.get("backbone", {}).values())
        matrix = P.panel_similarity_matrix(STATE["data"], labels)
        matrices.append({
            "ci": complexity if math.isfinite(complexity) else None,
            "fluors": [{
                "name": fl,
                "backbone": fl in backbone_fluors,
                "status": ("unmapped" if fl in panel["unmapped_fluors"]
                           else "substituted" if fl in panel["substitutions"] else ""),
                "note": panel["substitutions"].get(fl, ""),
            } for fl in labels],
            "values": [[round(float(v), 4) if math.isfinite(v) else None for v in row]
                       for row in matrix],
        })

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
                titles = []
                if ab in panel.get("backbone", {}):
                    titles.append("backbone")
                if fl in panel["unmapped_fluors"]:
                    cls = "unmapped"
                    titles.append("unmapped")
                elif fl in panel["substitutions"]:
                    sub = html_lib.escape(panel["substitutions"][fl])
                    cls = "substituted"
                    titles.append(f"substituted: {sub}")
                else:
                    cls = "assigned"
                title = f' title="{"; ".join(titles)}"' if titles else ""
                parts.append(f'<td class="{cls}"{title}>{esc(fl)}</td>')
            else:
                parts.append('<td class="unassigned" title="unassigned"></td>')
        parts.append("</tr>")

    parts.append('</tbody><tfoot><tr class="ci-row"><th colspan="2">Matrix_CI</th>')
    for c, _, _ in best:
        parts.append(f'<td class="num">{c:.4f}</td>')
    parts.append("</tr></tfoot></table>")

    substitutions = sorted({
        (fluor, substituted)
        for _, _, panel in best
        for fluor, substituted in panel["substitutions"].items()
    })
    substitution_parts = []
    if substitutions:
        substitution_parts.append('<div class="substitutions-heading">Substitutions made for calculations:</div>')
        for fluor, substituted in substitutions:
            substitution_parts.append(
                f'<div class="substitution-row">{esc(fluor)} &rarr; {esc(substituted)}</div>'
            )

    summary = (
        f"Maximum panel size: {STATE['max_count']}/{STATE['total']}. "
        f"Fluorophore sets scored: {evaluated:,} (branches pruned: {pruned:,})."
    )
    return json.dumps({
        "summary": summary,
        "html": "".join(parts),
        "substitutions_html": "".join(substitution_parts),
        "matrices_json": json.dumps(matrices, allow_nan=False),
    })


__export__ = ["setup", "search", "finalize"]
