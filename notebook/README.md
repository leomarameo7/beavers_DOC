# Notebook folder guide

This folder holds every Quarto document and R script behind the analysis. It is split into four
subfolders by what each file is *for*, not by date or topic:

| Folder | What's in it | Do these produce manuscript output? |
|---|---|---|
| [`manuscript/`](manuscript) | The documents that generate the figures, tables and numbers currently used in the manuscript (main text and Appendix S1). | **Yes — this is the current, load-bearing analysis.** |
| [`fitting/`](fitting) | R scripts that fit `brms` models, run cross-validation / hold-out tests, or run sensitivity analyses. No figures live here; each script caches its model or its numbers to `results/` and is read by a document in `manuscript/` or `correspondence/`. | Indirectly — these produce the cached results the documents above read. |
| [`correspondence/`](correspondence) | Analysis documents written to answer a specific co-author question or evaluate an option, not to be part of the manuscript itself. | No. |
| [`legacy/`](legacy) | The original, pre-revision reviewer notebook. Partly superseded, partly still current — see the note in that folder. | Partly — see below. |

Shared resources used by several documents stay at the top of `notebook/` rather than being copied
into each subfolder: `references.bib`, `ecology-letters.csl` (used by the documents in `manuscript/`
and `correspondence/`, referenced from one directory up in their YAML header), and `custom.css`
(currently unused by any document).

## `manuscript/` — current manuscript analysis

| Document | What it produces |
|---|---|
| `analysis_plan.qmd` | The main-text figures (1–3) and most of Appendix S1 (variable definitions, Q1–Q5, model validation, structural/causal sensitivity analysis, pathway decomposition, the methods map, the figure/table inventory). This is the single most up-to-date document and the one to read first. |
| `concentrations_fit.qmd` | The per-concentration measurement-error model that Figure 1b and Figure 3 of `analysis_plan.qmd` are built on. |
| `connectivity.qmd` | The alternative structural causal model with lateral stream–wetland connectivity as a driver of water residence time (Appendix S1, Figures S5–S6). |

## `fitting/` — model fitting, validation, sensitivity

Grouped by which extended or alternative model each script belongs to:

- **The two current SCM parameterizations** (`Q2_delta` / `Q2_downstream`, with measurement error on
  upstream DOC and the upstream → producer links): `me_refits.R` (fits both), `me_predictions.R`
  (leave-one-out cross-validation), `me_holdout.R` (80/20 stream hold-out test). Read by
  `correspondence/scm_outcome_choice.qmd`.
- **The plain downstream-DOC SCM** (no upstream measurement error, no producer links — the version
  `analysis_plan.qmd` currently reports): `scm_downstream_fit.R` (fit), `scm_downstream_fixed_fit.R`
  (robustness check: inheritance slope fixed at 1 instead of estimated), `scm_holdout.R` (80/20
  hold-out), `scm_predictive_checks.R` (leave-one-out).
- **Sensitivity and mechanism analyses** (section 9 and Q5 of `analysis_plan.qmd`):
  `dag_sensitivity.R` (alternative DAGs, prior sensitivity), `pathway_decomposition.R` (Q5 pathway
  decomposition of the local modification), `derived_ddoc_by_dag.R` (does the derived ΔDOC depend on
  the SCM specification — feeds `correspondence/scm_outcome_choice.qmd`).

None of these take long to fit from scratch except the SCM refits (tens of minutes); all cache their
fitted model or output to `results/models_fit/` or `results/*.rds` via `here()`, so re-running a
document that reads their output does not refit anything unless the cached file is deleted.

## `correspondence/` — co-author decision documents

- `scm_outcome_choice.qmd` — written for Francesco: compares the two SCM parameterizations
  (`Q2_delta`, `Q2_downstream`) side by side, including the extended (measurement-error) refits in
  `fitting/me_*.R`. Not part of the manuscript.
- `scm_rebuild.qmd` — the earlier decision document that first tested Josh's proposal (downstream DOC
  instead of ΔDOC as the SCM's outcome) on the real data with a cost/risk assessment; superseded by
  the fuller comparison in `scm_outcome_choice.qmd`, kept for the record.

## `legacy/` — the original reviewer notebook

`notebook.qmd` is the original, single-file notebook described in the repository's top-level
`README.md` ("the entire analysis"). It is **mixed status**, not simply superseded:

- Its **beaver amplification-factor (*M*) section and old Figure 2** are explicitly removed from the
  manuscript — `manuscript/analysis_plan.qmd` states this directly ("The amplification-factor section
  and its figure are removed; Q3 is answered by $R_s$ ... and $G$ ...").
- Its **Figures S1–S4, S7 and S8** (priors, posterior predictive checks, stream-level classification,
  total-effects path tracing) are listed as still-current, "existing" Appendix S1 figures in
  `analysis_plan.qmd`'s own figure inventory, and are not reproduced anywhere else.

So this file cannot simply be archived or deleted without losing figures the current Appendix S1
still relies on, but it also documents an analysis (*M*) that the manuscript no longer uses. It is
kept as-is, pending a decision on whether to split the still-current figures out into
`manuscript/`. Flag this to Leo before editing its content.
