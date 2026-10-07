# Notebook folder guide

Every Quarto document and R script behind the manuscript and the supporting information. The
top-level [`README.md`](../README.md) has the figure-by-figure table (which notebook draws which
figure and what it needs) and the order in which to run things.

| Folder | What's in it |
|---|---|
| [`manuscript/`](manuscript) | Quarto documents that draw the figures and compute the numbers of the manuscript and SI. Each figure is written directly to `overleaf/` by the chunk that draws it. |
| [`fitting/`](fitting) | R scripts that fit the heavy `brms` models and cache their results to `results/`. No figures are drawn here. |

Shared files at the top of `notebook/`: `references.bib` and `ecology-letters.csl` (bibliography and
citation style used by the notebooks, referenced from their YAML header).

## `manuscript/`

| Document | Draws | Reads |
|---|---|---|
| `analysis_plan.qmd` | Figures 1–4; SI Figures S7, S8, S11–S14; Table S2. Also the numbers for the hold-out, intensity and pathway-decomposition analyses | `data/processed/*`, `data/elevation/`, `results/me_refits.rds`, `results/dag_sensitivity_me.rds`, `results/me_predictions.rds`, `results/me_holdout.rds`, `results/pathway_decomposition.rds`, `results/concentrations_thin.rds` |
| `appendix_s1_figures.qmd` | SI Figures S1, S3–S6 (priors, residuals, posterior predictive check, posterior parameters, stream-level ΔDOC) | `results/models_fit/delta_SCM_me.rds` |
| `connectivity.qmd` | SI Figure S10 (lateral stream–wetland connectivity) | `m6.csv`, `later_connectivity_berger.csv` |
| `concentrations_fit.qmd` | The per-stream concentration posteriors (`results/concentrations*.rds`) used by Figure 4 | `m6.csv` |

## `fitting/`

| Script | Produces |
|---|---|
| `me_refits.R` | The two main SCMs (`delta_SCM_me`, `downstream_DOC_me`) and `results/me_refits.rds` |
| `me_predictions.R` | `results/me_predictions.rds`: predictions and leave-one-out for both models |
| `me_holdout.R` | `results/me_holdout.rds`: 80/20 hold-out by stream |
| `dag_sensitivity_me.R` | `results/dag_sensitivity_me.rds`: alternative DAGs and halved/doubled priors |
| `downstream_DOC_me_fixed_fit.R` | The fit with the inheritance slope fixed at 1 (robustness check in `analysis_plan.qmd`, section 4) |
| `scm_downstream_fit.R` | The plain downstream-DOC SCM, input of `pathway_decomposition.R` |
| `pathway_decomposition.R` | `results/pathway_decomposition.rds`: pathway decomposition of the local modification |

Fits use `seed = 7` and 4 chains. `results/models_fit/` is not tracked by git (size); every script
caches its model there, so rerunning a script or rendering a document reuses it.
