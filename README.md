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
- "Fit at most N antibodies" limit (WebUI, or `max_assigned` in Python): searches for the best panels with fewer antibodies than the maximum. This is no slower than the full search, usually faster
- WebUI: suggests the best 1 to 5 extra fluorophores to add to a chosen result panel (`suggest_additions` in Python), ranked by the resulting complexity, with the extended similarity matrix
- WebUI warning for fluorophore pairs with a similarity of 0.9 or more in the shown panel, with advice on how far to lower the limit
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

```bash
uv sync

# Save your input data in Data/input.csv, then
uv run .\Panel_optimized.py
```

The WebUI can be run locally at http://localhost:8000/: 

```bash
python -m http.server 
```
Requires Python 3.12+.
