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