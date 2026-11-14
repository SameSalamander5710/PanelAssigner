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

- Similarity matrices: a purple line separates the fluorophores you fixed
  (the optional input) from the ones the search chose, as a green line
  separates suggested additions. The purple line only appears when a panel
  has fixed fluorophores.
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
