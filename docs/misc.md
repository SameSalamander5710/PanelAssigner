## Miscellaneous

### Fluorophore aliases

Fluorophores may be written in alternative spellings. Before matching, names
are rewritten using the rules in `constants/Fluorophore_Aliases.csv`
(`Alias,Canonical`), in the CLI and in the web UI alike. For example
`Brilliant Violet 421` -> `BV421`, `AF647` -> `Alexa Fluor 647`,
`PB` -> `Pacific Blue`, `FVS450` -> `Fixable Viability Stain 450` and
`PE-AF700` -> `PE-Alexa Fluor 700`. Names are shown in their canonical form in
the results.

- **Canonical** is the common form that all spellings are rewritten to. The
  input, `Similarity_Index.csv` and `Detector_Fluorophores.csv` all go through
  the same rules, so the reference files may use any spelling (a matrix entry
  called `Brilliant Violet 421` matches an input of `BV421`). Two matrix rows
  that mean the same fluorophore give an error.
- A trailing `{n}` stands for the number (`421`, `5.5`, `575V`) and must appear
  on both sides: `AF {n},Alexa Fluor {n}` covers every Alexa Fluor.
- Matching ignores case, spaces, hyphens and slashes and works on whole tokens
  inside tandems. `(R)` / `(TM)` marks from vendor names are ignored.
- If an alias appears twice, the first row wins. Add your own rows at the end.

### Detector_Fluorophores.csv

Columns: `Detector,Fluorophore,Type`. `Detector` is the peak channel on the
Cytek Aurora 5L, and fluorophores sharing a channel can substitute for one
another when an input fluorophore is not in `Similarity_Index.csv`. When
several matrix fluorophores share a channel, the first one in file order is
used, so keep each channel's rows sorted.

`Type` says what the reagent is. Suggested additions (`suggest_additions`, and
the web UI) only come from rows whose Type is in `SUGGESTED_TYPES` in
`Panel_optimized.py`, by default `Antibody`:

- `Antibody`: sold as an antibody or streptavidin conjugate for flow cytometry
- `Dye`: viability, proliferation, cell-tracking and organelle stains
- `Fluorescent Protein`: genetically encoded reporters (GFP, mCherry, ...)
- `Microscopy`: label dyes that are mainly sold for microscopy

A fluorophore with no Type, or a file without a Type column, is not filtered.
A fluorophore you put in your own panel is always used, whatever its Type.

### Performance notes

- `prepare_reference()` applies the alias rules to the similarity matrix,
  detector map and Type table once. `worker.py` does this at start-up and
  passes `canonical=True` to `prepare_search`, `suggest_additions` and
  `eligible_fluorophore_count`; without the flag they do it on every call.
- `seed_bound()` (run once per worker in `setup`) keeps the similarity matrix and
  the antibody matching of the panel it is building and adds one fluorophore at a
  time, like the search does, instead of rebuilding both from names for every
  candidate. It returns the same bound as the plain version (`tests/test_seed_bound.py`
  compares them) in less than half the time on large panels.
- Compute and Suggest never run together: worker 0 holds the panels that Suggest reads
  and Compute rebuilds, so each button is disabled while the other runs (`active` in
  `main.py`).
- The panel search tries the most-overlapping fluorophores first
  (`search_order`), which cuts off more of the search tree. This changes how
  many sets are scored and pruned, not the panels found: every panel,
  complexity and tie-break uses the alphabetical `all_fluors` order.
- The metric scores (`_a_score` and friends) and `_score_matrix` run once per
  search node, tens of thousands of times for a panel of this size. They use
  `sum` and `math.*` on the eigenvalues instead of `np.mean`, `np.sqrt` and
  `np.exp`, and skip an `isfinite` pass (the matrix was already checked), which
  made the A- and D-optimality searches 15-30% faster. Most of what is left
  is `numpy.linalg.eigvalsh` itself.
- A set with two identical spectra has an infinite score (any metric) in
  theory; a non-positive or near-zero eigenvalue is reported as infinite. Such
  panels only appear when fewer than `top_n` ordinary panels exist.

### Input names and unmapped fluorophores

- Antibody names that differ only in case or spacing (`CD3` and `cd3`) are one
  antibody: the first spelling is kept and the fluorophore options of the
  others are added to it. This is the same rule the backbone already used.
- A fluorophore that is neither in `Similarity_Index.csv` nor shares a detector
  with a fluorophore that is, is *unmapped*. It is assumed to overlap with
  nothing, so panels containing it get a score that is better
  than the real one. Such cells are yellow and the results table has an `Unmapped` row with
  the count per panel. "Only use fluorophores that are in the matrix"
  (`require_mapped=True`) leaves them out instead; backbone fluorophores are
  never left out.

### Reading the page

- The page follows the order of the work. A step bar (sticky on wide screens)
  lists 1 Input, 2 Settings, 3 Results, 4 Inspect and 5 Extend, marks the one
  on screen, and says what to do next; Inspect and Extend stay locked until
  there are results. Input and Settings sit side by side on wide screens,
  Results use the full width, and Inspect and Extend sit side by side on
  screens of about 1360 px or more. New results are scrolled into view unless
  they are already visible. The selected "Rank by" score is described under the
  box, and the arrow buttons beside the panel and suggestion boxes step through
  them. These guides are plain JavaScript in the last inline script before the
  PyScript tags; `main.py` and the other scripts rely on the element ids, so
  keep them when changing the layout.
- Similarity matrices: a purple line separates the fluorophores you fixed
  (the optional input) from the ones the search chose, as a green line
  separates suggested additions. The purple line only appears when a panel
  has fixed fluorophores.
- Suggested additions are a table with a column per suggestion and a row per
  added fluorophore, then the score and its change from the panel. Clicking a
  column (or using the arrows and the number box) shows that suggestion's
  matrix, the same way the results table selects a panel.
- Leaving fluorophores out of the suggestions: the box above the Suggest
  button takes one name per line, or comma, semicolon or tab separated
  (`parse_exclusions`; `parseExclusions` in `index.html` does the same). They
  are only taken out of the pool of candidates, so the suggested sets are the
  best of what is left, as if those fluorophores were not in the reference
  matrix (`exclude=` of `suggest_additions`). Fluorophores of the panel itself
  are never affected. Names ignore case, spaces and slashes. `*` stands for any
  characters, anywhere: `BV*`, `PE-*`, `*605`. A pattern is compared with
  every spelling of a fluorophore: its canonical name and what each alias
  rule gives for it, so `BV*` and `Brilliant Violet*` both leave out BV421,
  and `AF700` leaves out Alexa Fluor 700. The flip side is that a pattern can
  reach a fluorophore through an unusual spelling (`R*` also takes the "R-PE-"
  forms of the PE tandems). After Suggest, the note under the box names what
  each pattern removed (hover for the full list) and warns about a pattern
  that matched nothing, which is usually a typo.
- A minus appears on a suggested fluorophore when the pointer is on it (when it
  has keyboard focus, and always on touch screens). Clicking it adds the
  name to the box; it is struck through in the table, and a bar offers "Suggest
  again", because the suggestions on screen were made without it. Editing the
  box has the same effect. The box keeps its text when the panels are
  computed again or another panel is selected.
- Several assignments per panel: the score depends only on the set of
  fluorophores, so antibodies that share two or more fluorophores can swap
  them without changing it (`CD45` and `CD11b`, both `FITC` or `BV421`, give
  two assignments of the same panel). `panel_assignments(data, panel)` lists
  every way to give the panel's fluorophores to distinct antibodies. With
  "fit at most", which antibodies are left out varies as well; the fixed
  backbone never does. The first is the one the search reports
  (`panel["assignments"]`), the others follow in input order (the first
  antibody's first-listed fluorophore first). At most `MAX_ALTERNATIVES` (200)
  are kept per panel, because the number can explode (6 antibodies that share
  6 fluorophores have 720); the count is then shown as `200+`. The listing
  keeps a matching of the fluorophores still to place and only follows moves
  that can be completed, so it takes time per assignment found, not per
  possible assignment (a few milliseconds for 200 in CPython, on panels of up
  to 25 antibodies). It runs in `finalize` for the shown panels only, not in
  the search. `find_best_panels` does not call it: its `ranked_panels` keep
  the one assignment they always had.
  A panel with more than one gets arrows and a counter under its name in the
  results table; ordinary results look as before. Stepping (it wraps around)
  rewrites that column of the table and selects the panel, and the line
  above the matrix names the assignment (`A-optimality: 1.0000 · assignment 2
  of 24`) and the antibody labels in its warnings follow it. The matrix itself
  keeps its order, so it does not reshuffle. Each assignment goes to the page
  as a list of table rows, one per fluorophore of the panel, in the
  `assignments` of `matrices_json` (`antibodies` names the rows).
  `_result_cell` (`worker.py`) and `paintAssignment` (`index.html`) must give
  the same cell. The arrows are marked `nocopy`, and Copy takes the table as
  shown, so it copies the assignment on screen.
- Clicking a panel column in the results table selects that panel, the same
  as typing its number in the matrix panel box. The selected column is marked
  with a style sheet, so the table that Copy takes is unchanged. Clicks that
  end a text selection are ignored.
- Above Compute, a live summary (`worker.summarize`) checks the inputs as you
  type, with no search: antibodies, pairs, fixed fluorophores, and how many
  fluorophores are detector-substituted (blue) or unmapped (yellow, "left
  out" if only fluorophores in the matrix are used). The same colours mark
  the rows of both text boxes. Each text box has a copy of its text behind it
  in transparent ink so that rows can be coloured, so font, line height,
  padding and wrapping of `.hl-back` must stay identical to those of the
  `textarea` (see the CSS comment).
  The check runs in `refresh_live` (`main.py`) as a single loop (`live_loop`,
  a task) with at most one call to its worker at a time: editing while a call
  is running does not start another one, the loop checks again with the new
  text when the call returns. A check that takes longer than `LIVE_GIVE_UP`
  seconds is abandoned, and an error from the worker is shown in the summary.
  Compute never waits for a check. It starts at once on the workers that
  have no call out (`busy`, kept by `call`), and a worker that is still
  answering a check joins in as soon as it has answered, or is left out if
  the others finish first or it never answers. The answer of a check that
  was already with a worker when Compute was clicked is dropped, and the
  summary is refreshed when Compute is done (`refresh_live(True)`). The
  check uses the last free worker, never the first (it holds the panels
  for Suggest).
- **Nothing that JavaScript starts may suspend (`main.py`).** MicroPython in
  the browser runs a coroutine that JavaScript starts (an `async` `@when`
  handler, a top-level `await`) as a "top-level" coroutine, and all of those
  share one wake-up slot (`TopLevelCoro` in MicroPython's
  `ports/webassembly/asyncio/core.py`). When one is suspended in a sleep, a
  worker call or `gather` and another starts and suspends, they overwrite each
  other: one is never resumed (a live summary that never arrives, a Compute
  that waits for 15 s), the scheduler throws into the shared slot
  (`AttributeError: type object 'TopLevelCoro' has no attribute 'throw'`), or
  the CPU spins (typing into the boxes froze the page). So event handlers are
  plain functions that do their synchronous work and hand the rest to
  `spawn()`, which makes a task with an identity of its own, and the module
  ends with `spawn(start())` instead of `await start()`. CPython's asyncio
  does not behave like this, so `tests/test_live_check.py` enforces the rule
  directly (`entry_point_problems`: no `async` `@when` handler, no top-level
  `await`) and calls handlers the way JavaScript does.
- "Load example panel" fills both text boxes from `constants/Example_Panel.csv`
  and `constants/Example_Backbone.csv` (12 human immunophenotyping markers
  with 2-3 fluorophore options each, and an `FVS 700` viability dye to avoid).
  It is plain JavaScript in `index.html`, so it works while the Python workers
  are still loading; the live summary appears once they are ready. Several
  names are written as aliases on purpose (`BUV 395`, `Alexa 700` and `AF 700`
  are the same dye, `PB`, `R-PE`, ...). The files are fetched from the site
  itself, so the button needs the page to be served over http(s), not opened
  from disk. If you edit them, `tests/test_example_panel.py` checks that they
  still read cleanly (nothing substituted or unmapped) and give a full panel.
- "Fit at most N of M antibodies" counts only the antibodies of the main
  text box (M is the live summary's antibody count). The fixed ones are always
  in the panel and are not counted, so a panel in the results table has up to
  N plus the number of fixed rows. In Python this is `max_free`;
  `max_assigned` is the whole panel size, fixed rows included. The "lower
  Fit at most to ..." advice under a similarity matrix uses the same counting
  and never suggests leaving out a fixed fluorophore.
