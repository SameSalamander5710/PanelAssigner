# PanelAssigner

Simple optimizer for flow cytometry panel design. From a pool of antibody-fluorophore pairs, optimally assigns fluorophores to antibodies and ranks panel configurations by spectral complexity. 

Mainly intended to optimize panel design from a constrained list of antibodies. Does NOT replace optimization based on biology and gating strategies.

The tool is hosted online at https://samesalamander5710.github.io/PanelAssigner/

## Features

- Evaluates all feasible fluorophore combinations using exact search with pruning
- Ranks panels by spectral overlap, using a metric of your choice: A-, D- or E-optimality, or the condition index of the similarity matrix (dropdown in the WebUI, `metric=` in Python)
- "Load example panel" button (WebUI) that fills both text boxes with a 12-marker immunophenotyping example, so you can try the tool straight away
- Add fluorophores to penalize overlap with (for when you want to add onto an existing panel, avoid autofluorescence, or stain labelled/tagged cells)
- Detector-based substitution for fluorophores not in the reference matrix
- Similarity matrix plot of any result panel (WebUI only), below the results
- Suggests the best 1-5 extra fluorophores to add to a chosen result panel (WebUI; `suggest_additions` in Python), ranked by the resulting complexity, shown as a table with a column per suggestion
- Fluorophores can be left out of the suggestions (not for sale, or just not wanted): type them in the WebUI, with `*` as a wildcard (`BV*`), or click the minus that appears on a suggested fluorophore (`exclude=` in Python)

## Input Format

CSV or text pasted from Excel, with two columns:

| Antibody | Fluorophore |
|----------|-------------|
| CD3      | BV421       |
| CD3      | FITC        |
| CD4      | APC         |

The WebUI's "Load example panel" button uses `constants/Example_Panel.csv` and `constants/Example_Backbone.csv` in this format; they are also a template for your own input.

## Output

Each result includes:
- **Assignments**: antibody-to-fluorophore mappings
- **Score** (A-optimality, D-optimality, E-optimality or condition index, as chosen): lower values indicate less spectral overlap
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

## Complexity calculations

Panels are ranked by one of four scores of the panel's similarity matrix (`METRICS` in `Panel_optimized.py`; `metric=` argument of `find_best_panels` and `suggest_additions`; "Rank by" dropdown in the WebUI). The default is `DEFAULT_METRIC` (A-optimality).

For a panel of $n$ fluorophores, build the symmetric matrix $S$ ($n \times n$) where $S_{ii}=1$, $S_{ij}$ is the spectral similarity between fluorophores $i$ and $j$ from `constants/Similarity_Index.csv`, and cells involving an unmapped fluorophore are $0$. Taking the eigenvalues $\lambda_i$ of $S$:

| Metric (`metric=`) | Score | Meaning |
| --- | --- | --- |
| Condition index (`condition_index`) | $\sqrt{\lambda_{\max}/\lambda_{\min}}$ | Numerical stability of the matrix; ignores everything but the two extreme eigenvalues |
| A-optimality (`a_optimality`) | $\frac{1}{n}\sum_i \frac{1}{\lambda_i}$ | Average variance inflation (mean VIF): every overlapping fluorophore counts |
| D-optimality (`d_optimality`) | $\det(S)^{-1/n}=\left(\prod_i \frac{1}{\lambda_i}\right)^{1/n}$ | Total information: the volume spanned by the spectra (geometric mean of $1/\lambda_i$) |
| E-optimality (`e_optimality`) | $1/\lambda_{\min}$ | Worst-direction inflation: only the most overlapping combination counts |

Lower values are better for all four, and a perfectly orthogonal panel of any size scores $1.0$. If the smallest eigenvalue is $\le 10^{-12}$ or the matrix contains non-finite values, the score is reported as infinite. Panels with fewer than two fluorophores score $1.0$. For any panel, E $\ge$ A $\ge$ D $\ge 1$ (maximum, arithmetic mean and geometric mean of the $1/\lambda_i$).

The A-, D- and E-optimality criteria follow Rajwa & Roederer, *Lost in Translation: Harmonizing Terminology and Defining Mathematical Tools for Panel Optimization*, Cytometry A 107(12), 2025 ([doi:10.1002/cyto.a.24972](https://doi.org/10.1002/cyto.a.24972)). There, the cosine-similarity matrix of $\ell_2$-normalized spectra acts as a Fisher-information-like matrix, its inverse bounds the unmixing variances (Cramér-Rao), A-optimality is the trace of that inverse, D-optimality its determinant, and E-optimality its largest eigenvalue. Two differences here: the scores are divided by or rooted to $n$ so that panels of different sizes (backbone included) can be compared, which does not change the ranking among panels of the same size, and only the similarity index is used, not a full mixing matrix or a noise model. D-optimality is reported as $\Delta^{-1}$ of the paper's scaled measure $\Delta=\det(S)^{1/n}$ so that lower is better, like the others.

The search stays exact for every metric: each score has a proven lower bound for the final panel that can be computed from a partial set (see `_Metric` in `Panel_optimized.py`), so branches that cannot reach the top results are skipped.
