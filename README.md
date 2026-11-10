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
- Optional fixed backbone panel that is always included and takes priority in clashes
- Parallel execution across multiple CPU cores

## Input Format

CSV or pasted text with two columns:

| Antibody | Fluorophore |
|----------|-------------|
| CD3      | BV421       |
| CD3      | FITC        |
| CD4      | APC         |

### Fixed backbone (optional)

A backbone panel can be fixed in place with the same two-column format. Its
assignments are always included, are used in the complexity calculation, and
take priority when they clash with an input antibody or fluorophore.

- CLI: save the pairs in `Data/input_backbone.csv` (header optional).
- Web UI: paste them into the **Fixed Panel Backbone** box.

## Output

Each result includes:
- **Assignments**: antibody-to-fluorophore mapping
- **Unassigned antibodies**: those that could not be matched
- **Matrix Complexity**: lower values indicate less spectral overlap
- **Substitutions**: fluorophores replaced by detector-matched alternatives
- **Unmapped**: fluorophores not found in the reference matrix

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
