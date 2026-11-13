# PanelAssigner

Simple optimizer for flow cytometry panel design. Assigns fluorophores to antibodies and ranks panel configurations by spectral complexity. 

Mainly intented to optimize panel design from a constrained list of antibodies. Does NOT replace optimization based on biology and gating strategies.

The tool is hosted online at https://samesalamander5710.github.io/PanelAssigner/

## Features

- Finds the maximum number of antibodies that can be assigned distinct fluorophores
- Evaluates all feasible fluorophore combinations using exact search with pruning
- Ranks panels by spectral overlap (A-optimality of the similarity matrix)
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
- **A-optimality**: lower values indicate less spectral overlap
- **Substitutions**: fluorophores replaced by detector-matched alternatives, if not found in the current database
- **Unmapped**: fluorophores not found in the reference matrix
- **Similarity matrix**: similarity matrix of result panels

## Complexity calculations

The A-optimality score (implemented in `_complexity_for_set`) is the mean variance inflation factor of the panel's similarity matrix.

For a panel of $n$ fluorophores, build the symmetric matrix $S$ ($n \times n$) where $S_{ii}=1$, $S_{ij}$ is the spectral similarity between fluorophores $i$ and $j$ from `constants/Similarity_Index.csv`, and cells involving an unmapped fluorophore are $0$. Taking the eigenvalues $\lambda_i$ of $S$:

$$
\boxed{\text{Score}_{A}(S)=\frac{1}{n}\sum_{i=1}^{n}\frac{1}{\lambda_i}}
$$

Lower values are better: a perfectly orthogonal panel of any size scores $1.0$, and the score grows as overlaps make the matrix closer to singular. If any eigenvalue is $\le 0$ (or below a numerical tolerance of $10^{-12}$) or the matrix contains non-finite values, the score is reported as infinite. Panels with fewer than two fluorophores score $1.0$.


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
