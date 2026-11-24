"""One of several identical PyScript (Pyodide) workers. Each runs in its own
Web Worker thread; main.py hands out search subtrees to whichever is free."""
import csv
import html as html_lib
import json
import math

import Panel_optimized as P

SIM_DF, SIM_NAMES = P.load_similarity_matrix("./constants/Similarity_Index.csv")
DET_MAP = P.load_detector_map("./constants/Detector_Fluorophores.csv")
TYPES = P.load_fluorophore_types("./constants/Detector_Fluorophores.csv")  # only Antibody rows are suggested
try:
    ALIASES, ALIAS_ERROR = P.load_fluorophore_aliases("./constants/Fluorophore_Aliases.csv"), None
except Exception as exc:  # a typo in the CSV is reported in the UI, not lost in the worker
    ALIASES, ALIAS_ERROR = [], str(exc)

# Apply the alias rules to the reference data once, not on every Compute and Suggest.
try:
    REF_SIM, REF_NAMES, REF_DET, REF_TYPES = P.prepare_reference(
        SIM_DF, SIM_NAMES, DET_MAP, ALIASES, TYPES
    )
    REF_ERROR = None
except ValueError as exc:  # e.g. two matrix names that mean the same fluorophore
    REF_SIM = REF_NAMES = REF_DET = REF_TYPES = None
    REF_ERROR = str(exc)

TARGET_TASKS = 32  # more, smaller subtrees = better load balancing and pruning
STATE = {}


def _line_delimiter(line):
    """Separator of one pasted row: tab, comma, else semicolon as fallback.

    A separator inside a quoted field is part of the value, not a separator,
    so a comma in a quoted name does not turn an Excel row into a CSV row.
    """
    quoted = False
    has_comma = False
    for char in line:
        if char == '"':
            quoted = not quoted
        elif not quoted:
            if char == "\t":
                return "\t"
            if char == ",":
                has_comma = True
    return "," if has_comma else ";"


def _pasted_rows(text):
    """Yield (line number, fields) from a pasted text, with a separator per row.

    Typed CSV rows and rows pasted from Excel can then share one text box.
    Line numbers count the lines as they stand in the box, including any
    stripped off the top by .strip().
    """
    body = text.strip()
    skipped = text[: len(text) - len(text.lstrip())].count("\n")
    for index, line in enumerate(body.split("\n")):
        # Trailing whitespace (a stray tab, spaces, the \r of CRLF) must not
        # decide the separator of a row such as "CD3,BV421<tab>".
        line = line.rstrip()
        parsed = next(csv.reader([line], delimiter=_line_delimiter(line)), [])
        yield skipped + index, parsed


def parse_pasted(text, required=True):
    """Parse Antibody/Fluorophore pairs (tab, comma or semicolon, per row)."""
    antibodies = {}
    for _, row in _pasted_rows(text):
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
    """Parse optional fixed backbone entries; the first row per antibody wins.

    A row with a single value is a bare fluorophore and gets an automatic
    antibody name (Avoid_1, Avoid_2, ...); rows with two values are treated as
    Antibody/Fluorophore pairs."""
    text = text or ""
    if not text.strip():
        return {}
    backbone = {}
    auto = 0
    for _, row in _pasted_rows(text):
        fields = [field.strip() for field in row if field.strip()]
        if not fields:
            continue
        if len(fields) >= 2:
            ab, fl = fields[0], fields[1]
            if ab.lower() == "antibody" and fl.lower() == "fluorophore":
                continue
        else:
            auto += 1
            ab, fl = f"Avoid_{auto}", fields[0]
        if ab not in backbone:
            backbone[ab] = fl
    if not backbone:
        raise ValueError("No Antibody/Fluorophore pairs found in the fixed backbone.")
    return backbone


def _row_fluors(text, backbone=False):
    """[(line number, fluorophore)] for each usable row, read line by line.

    This follows parse_pasted (backbone=False) and parse_backbone (True), and
    numbers the lines as they stand in the text box, so rows can be highlighted."""
    rows = []
    for number, parsed in _pasted_rows(text or ""):
        if not parsed:
            continue
        if backbone:
            fields = [field.strip() for field in parsed if field.strip()]
            if not fields:
                continue
            if len(fields) >= 2 and fields[0].lower() == "antibody" and fields[1].lower() == "fluorophore":
                continue
            rows.append((number, fields[1] if len(fields) >= 2 else fields[0]))
        else:
            if len(parsed) < 2:
                continue
            ab, fl = parsed[0].strip(), parsed[1].strip()
            if not ab or not fl or (ab.lower() == "antibody" and fl.lower() == "fluorophore"):
                continue
            rows.append((number, fl))
    return rows


def summarize(text, backbone_text="", require_mapped=""):
    """A quick look at the input, for showing before anything is computed.

    Counts antibodies and pairs as the search will see them, and flags
    fluorophores that are detector-substituted or unmapped, overall and per
    line of each text box (so the rows can be coloured). It uses the same
    alias, substitution and parsing code as Compute, but no search."""
    try:
        if ALIAS_ERROR:
            raise ValueError(f"Could not load fluorophore aliases: {ALIAS_ERROR}")
        if REF_ERROR:
            raise ValueError(REF_ERROR)
        require_mapped = str(require_mapped).strip().lower() in ("1", "true", "yes", "on")
        lines = {"input": text.split("\n"), "backbone": (backbone_text or "").split("\n")}
        rows = {"input": _row_fluors(text), "backbone": _row_fluors(backbone_text or "", True)}

        # Each distinct fluorophore (by canonical name) is looked up once.
        canonical = {}
        for kind in rows:
            for _, fluor in rows[kind]:
                canonical.setdefault(fluor, P.canonicalize_fluorophore(fluor, ALIASES))
        _, _, substituted, unmapped = P._prepare_similarity(
            REF_SIM, REF_NAMES, canonical.values(), REF_DET
        )
        def status(fluor):
            name = canonical[fluor]
            return "unmapped" if name in unmapped else ("substituted" if name in substituted else "")

        flagged = {"substituted": {}, "unmapped": {}}
        row_status = {kind: [""] * len(lines[kind]) for kind in rows}
        for kind in rows:
            for number, fluor in rows[kind]:
                kind_of = status(fluor)
                row_status[kind][number] = kind_of
                if kind_of:
                    note = substituted.get(canonical[fluor], "")
                    flagged[kind_of].setdefault(canonical[fluor], {"name": fluor, "note": note})

        try:
            main = parse_pasted(text, required=False)
            fixed = parse_backbone(backbone_text) if (backbone_text or "").strip() else {}
        except ValueError:
            main, fixed = {}, {}
        counts = None
        if main:
            try:
                data = P.prepare_search(
                    main, REF_SIM, REF_NAMES, REF_DET, 1, fixed, ALIASES,
                    canonical=True, require_mapped=require_mapped,
                )
            except ValueError:
                data = None
            if data is not None:
                counts = (
                    len(data["free_original_order"]),
                    sum(len(v) for v in data["antibodies_by_fluor"].values()),
                    len(data["backbone_assignments"]),
                )
        if counts is None:  # nothing assignable (yet): count what was typed
            counts = (len(main), sum(len(v) for v in main.values()), len(fixed))
        return json.dumps({
            "ok": True,
            "antibodies": counts[0], "pairs": counts[1], "fixed": counts[2],
            "require_mapped": require_mapped,
            "substituted": list(flagged["substituted"].values()),
            "unmapped": list(flagged["unmapped"].values()),
            "input_rows": row_status["input"],
            "backbone_rows": row_status["backbone"],
        })
    except Exception as exc:  # shown under the inputs, never lost in the worker
        return json.dumps({"ok": False, "error": str(exc)})


def setup(text, top_n=10, backbone_text="", max_fit="", require_mapped="",
          metric=P.DEFAULT_METRIC):
    """Every worker runs this with the same input and builds the same task list.

    `metric` is the key of the score panels are ranked by (`P.METRICS`).

    `max_fit` is the "fit at most N antibodies" limit on the antibodies of the
    main text box (the fixed backbone is always included); blank means no limit.
    `require_mapped` ("1", "true" or True) leaves out fluorophores that have
    no similarity data instead of treating them as overlapping with nothing."""
    try:
        if ALIAS_ERROR:
            raise ValueError(f"Could not load fluorophore aliases: {ALIAS_ERROR}")
        if REF_ERROR:
            raise ValueError(REF_ERROR)
        top_n = int(top_n)
        metric = str(metric or P.DEFAULT_METRIC)
        P.get_metric(metric)  # an unknown name is reported in the UI
        max_fit = str(max_fit).strip()
        if max_fit and not max_fit.isdigit():
            raise ValueError('"Fit at most" must be a whole number of antibodies.')
        max_free = int(max_fit) if max_fit else None
        require_mapped = str(require_mapped).strip().lower() in ("1", "true", "yes", "on")
        backbone = parse_backbone(backbone_text)
        antibodies = parse_pasted(text, required=not backbone)
        data = P.prepare_search(
            antibodies, REF_SIM, REF_NAMES, REF_DET, top_n,
            backbone=backbone, aliases=ALIASES, max_free=max_free,
            canonical=True, require_mapped=require_mapped, metric=metric,
        )
        if data is None:
            raise ValueError(
                "None of the fluorophores have similarity data, so nothing can be "
                "assigned. Untick the option to use only fluorophores in the matrix."
                if require_mapped else "No assignable antibody/fluorophore pairs."
            )
        tasks = P.make_search_tasks(data, TARGET_TASKS)
        STATE.update(data=data, tasks=tasks, top_n=top_n, metric=metric,
                     total=len(data["original_order"]),
                     free_total=len(data["free_original_order"]),
                     max_count=data["max_count"],
                     max_possible=data["max_possible"],
                     merged=data["merged_antibodies"],
                     dropped_unmapped=data["dropped_unmapped"])
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


def _result_cell(panel, antibody):
    """The table cell of `antibody` in `panel` (its fluorophore, or unassigned).

    The page repaints these cells when it steps to another assignment
    (`paintAssignment` in index.html); the two must give the same cell.
    """
    if antibody not in panel["assignments"]:
        return '<td class="unassigned" title="unassigned"></td>'
    fl = panel["assignments"][antibody]
    titles = []
    if antibody in panel.get("backbone", {}):
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
    return f'<td class="{cls}"{title}>{html_lib.escape(str(fl), quote=False)}</td>'


def _assignment_stepper(panel_no, count, more):
    """The arrows under a panel's name that step through its assignments.

    Empty for a panel with a single assignment, so ordinary results look as they
    always did. The page (index.html) wires the buttons up and keeps the
    counter; `nocopy` keeps this out of what Copy takes. `more` means that the
    search for assignments stopped at a limit, so the count is a minimum.
    """
    if count < 2:
        return ""
    plus = "+" if more else ""
    ways = f"at least {count}" if more else str(count)
    return (
        f'<span class="alt-nav nocopy" data-widest="{count}/{count}{plus}" '
        f'title="Same fluorophores, same score: {ways} ways to give them to the antibodies. '
        'Step through them.">'
        f'<button type="button" class="alt-btn" data-step="-1" '
        f'aria-label="Previous assignment of panel {panel_no}">&#9664;</button>'
        f'<span class="alt-count" aria-live="polite">1/{count}{plus}</span>'
        f'<button type="button" class="alt-btn" data-step="1" '
        f'aria-label="Next assignment of panel {panel_no}">&#9654;</button>'
        "</span>"
    )


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
    #
    # A panel is ranked by its fluorophores; when antibodies share fluorophores
    # the same set can be dealt out to them in several ways, all with this
    # score. Each way goes out as a list that gives, for each fluorophore of
    # the panel (same order as "fluors"), the index of its antibody in
    # "antibodies"; the first is the one the table shows. The matrix keeps the
    # order of the first, so it does not reshuffle when the page steps through.
    antibody_index = {ab: i for i, ab in enumerate(STATE["data"]["original_order"])}
    matrices = []
    for complexity, _, panel in best:
        labels = list(panel["assignments"].values())  # same order as the table
        antibody_of = {fl: ab for ab, fl in panel["assignments"].items()}
        backbone_fluors = set(panel.get("backbone", {}).values())
        matrix = P.panel_similarity_matrix(STATE["data"], labels)
        alternatives, more = P.panel_assignments(STATE["data"], panel)
        assignments = []
        for alternative in alternatives:
            holder = {fl: ab for ab, fl in alternative.items()}
            assignments.append([antibody_index[holder[fl]] for fl in labels])
        matrices.append({
            "score": complexity if math.isfinite(complexity) else None,
            "assignments": assignments,
            "more_assignments": more,
            "candidates": P.eligible_fluorophore_count(
                REF_SIM, REF_NAMES, REF_DET, labels, ALIASES, REF_TYPES,
                canonical=True,
            ),
            "fluors": [{
                "name": fl,
                "antibody": antibody_of[fl],
                "backbone": fl in backbone_fluors,
                "status": ("unmapped" if fl in panel["unmapped_fluors"]
                           else "substituted" if fl in panel["substitutions"] else ""),
                "note": panel["substitutions"].get(fl, ""),
            } for fl in labels],
            "values": [[round(float(v), 4) if math.isfinite(v) else None for v in row]
                       for row in matrix],
        })

    STATE["panel_labels"] = [[f["name"] for f in m["fluors"]] for m in matrices]

    parts = [
        '<table class="results-table"><thead><tr>',
        '<th class="num">#</th><th>Antibody</th>',
    ]
    for i, shown in enumerate(matrices, 1):
        parts.append(f"<th>Panel {i}{_assignment_stepper(i, len(shown['assignments']), shown['more_assignments'])}</th>")
    parts.append("</tr></thead><tbody>")

    for row_no, ab in enumerate(STATE["data"]["original_order"], 1):
        parts.append(f'<tr><td class="num">{row_no}</td><td>{esc(ab)}</td>')
        for _, _, panel in best:
            parts.append(_result_cell(panel, ab))
        parts.append("</tr>")

    label = P.get_metric(STATE["metric"]).label
    parts.append(f'</tbody><tfoot><tr class="ci-row"><th colspan="2">{esc(label)}</th>')
    for c, _, _ in best:
        parts.append(f'<td class="num">{c:.4f}</td>')
    parts.append("</tr></tfoot></table>")
    # Fluorophores without similarity data count as overlapping with nothing,
    # which flatters the score above; the hover notes say how many.
    n_unmapped = [len(panel["unmapped_fluors"]) for _, _, panel in best]

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

    if STATE["max_count"] < STATE["max_possible"]:
        size = (f"Panel size limited to {STATE['max_count']}/{STATE['total']} "
                f"({STATE['max_possible']} possible). ")
    else:
        size = f"Maximum panel size: {STATE['max_count']}/{STATE['total']}. "
    summary = (
        size +
        f"Fluorophore sets scored: {evaluated:,} (branches pruned: {pruned:,})."
    )
    notes = []
    with_unmapped = sum(1 for n in n_unmapped if n)
    if with_unmapped:
        notes.append(
            f"{with_unmapped} of {len(best)} panels contain fluorophores with no similarity "
            "data (yellow); their overlap is not counted, so their score is too good. "
            "Tick the option to use only fluorophores in the matrix to leave them out."
        )
    if STATE["dropped_unmapped"]:
        notes.append("Left out (no similarity data): " + ", ".join(STATE["dropped_unmapped"]) + ".")
    if STATE["merged"]:
        notes.append("Treated as the same antibody: " + ", ".join(
            f"{gone} = {kept}" for gone, kept in STATE["merged"]) + ".")
    return json.dumps({
        "summary": summary,
        "notes": notes,
        "has_unmapped": with_unmapped > 0,
        "html": "".join(parts),
        "substitutions_html": "".join(substitution_parts),
        "matrices_json": json.dumps({"total": STATE["total"],
                                     "free_total": STATE["free_total"], "panels": matrices,
                                     "antibodies": STATE["data"]["original_order"],
                                     "metric": STATE["metric"], "metric_label": label},
                                    allow_nan=False),
    })


def suggest(panel_no, n=1, exclude=""):
    """Best `n` fluorophores to add to result panel `panel_no` (1-based).

    `exclude` is the text of the "leave out" box: fluorophores that must not be
    suggested, one per line or comma separated, `*` as a wildcard (see
    Panel_optimized.parse_exclusions). Needs finalize() to have run on this
    worker. The names go back as JSON and the page inserts them as text.
    """
    try:
        patterns = P.parse_exclusions(exclude)
        panels = STATE.get("panel_labels")
        if not panels:
            raise ValueError("Compute panels first.")
        panel_no, n = int(panel_no), int(n)
        if not 1 <= panel_no <= len(panels):
            raise ValueError(f"There is no panel {panel_no}.")
        result = P.suggest_additions(
            REF_SIM, REF_NAMES, REF_DET, panels[panel_no - 1], n, top_n=5,
            aliases=ALIASES, types=REF_TYPES, canonical=True,
            metric=STATE["metric"],  # the metric the panels were ranked by
            exclude=patterns,
        )
        return json.dumps({
            "ok": True, "panel": panel_no, "n": n,
            "exclude": patterns,  # as applied; the page compares it with the box
            "excluded": result["excluded"],  # [{"pattern", "removed"}] for each
            "metric_label": P.get_metric(STATE["metric"]).label,
            "base_score": result["base_complexity"],
            "candidates": result["candidates"],
            "combinations": result["combinations"],
            "suggestions": [{"fluors": s["fluors"], "score": s["matrix_complexity"],
                             "values": s["values"]} for s in result["suggestions"]],
        }, allow_nan=False)
    except Exception as exc:
        return json.dumps({"ok": False, "error": str(exc)})


__export__ = ["setup", "search", "finalize", "suggest", "summarize"]
