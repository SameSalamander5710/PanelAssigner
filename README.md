# PanelAssigner

Simple optimizer for flow cytometry panel design. Assigns fluorophores to antibodies and ranks panel configurations by spectral complexity. 

Mainly intented to optimize panel design from a constrained list of antibodies. Does NOT replace optimization based on biology and gating strategies.

The tool is hosted online at https://samesalamander5710.github.io/PanelAssigner/

## Features

- Finds the maximum number of antibodies that can be assigned distinct fluorophores
- Evaluates all feasible fluorophore combinations using exact search with pruning
- Ranks panels by spectral overlap (condition index of the similarity matrix)
- Automatic fluorophore name matching against a reference similarity matrix
- Detector-based substitution for fluorophores not in the reference matrix
- Similarity matrix plot of any result panel (WebUI only), below the results
- Suggests the best 1-5 extra fluorophores to add to a chosen result panel (WebUI; `suggest_additions` in Python), ranked by the resulting complexity, with the extended similarity matrix. Only antibody-label fluorophores are suggested (not dyes, fluorescent proteins or microscopy labels), see `docs/misc.md`
- Parallel execution across multiple CPU cores

## Input Format

CSV or pasted text with two columns:

| Antibody | Fluorophore |
|----------|-------------|
| CD3      | BV421       |
| CD3      | FITC        |
| CD4      | APC         |

## Output

Each result includes:
- **Assignments**: antibody-to-fluorophore mapping
- **Unassigned antibodies**: those that could not be matched
- **Matrix Complexity**: lower values indicate less spectral overlap
- **Substitutions**: fluorophores replaced by detector-matched alternatives, if not found in the current database
- **Unmapped**: fluorophores not found in the reference matrix
- **Similarity matrix**: similarity matrix of result panels

## Installation

Clone the repo, then:

Run the WebUI locally at http://localhost:8000/: 

```bash
python -m http.server 
```

Or, to use the python scripts directly,

```bash
uv sync

# Save your input data in Data (input.csv, input_backbone.csv), then
uv run .\Panel_optimized.py
```

Requires Python 3.12+.
