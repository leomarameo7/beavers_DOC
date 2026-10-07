# Overleaf project: manuscript + supporting information

Import this GitHub repository into Overleaf (New Project → Import from GitHub).
Everything in the repo appears in the Overleaf file tree; only this folder is LaTeX.

| File | What it is |
|---|---|
| `main.tex` | The manuscript: Introduction, Results and Discussion, Methods. **Set as Main document in Overleaf.** |
| `references.bib` | Bibliography for the manuscript (`\cite{}` keys = first author + year, e.g. `larsen2021`) |
| `sn-jnl.cls`, `sn-nature.bst`, `sn-basic.bst` | LaTeX class and bibliography-style files needed to compile `main.tex` |
| `figures/media/` | Main-text figures `figure1_v2.png` (Fig. 1), `Figure_2a.png` (Fig. 2), `figure4_v1.png` (Fig. 3), `figure3_v2.png` (Fig. 4). Written by `notebook/manuscript/analysis_plan.qmd`; do not edit by hand |
| `supplementary/supporting_information.tex` | The Supporting Information, a standalone document. To compile it in Overleaf: Menu → *Main document* → choose this file |
| `supplementary/figures/media/` | SI figures, written by `analysis_plan.qmd`, `appendix_s1_figures.qmd` and `connectivity.qmd`, except `figS1_dag_hypotheses.png` and `figS8_dag_connectivity_alt.png`, which are drawn by hand |

The top-level `README.md` lists which notebook draws each figure.
