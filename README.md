# Pathways and magnitude of beaver-mediated DOC change in streams

Code and data for the manuscript **"Beaver engineering creates intense, context-dependent dissolved organic carbon control points in streams"**
(Capitani L. *et al.*).
Contact: Leonardo Capitani (leocapi07@gmail.com), WSL & Eawag, Switzerland.

## Research questions

Beaver dams turn flowing streams into ponds. We quantify how this beaver engineering changes
**dissolved organic carbon (DOC)** in Swiss streams sampled upstream and downstream of beaver
ponds in winter and summer (180 streams; 158 with both seasons complete are used for the
seasonal analyses), and ask:

1. **Direction** – Is the change in DOC across a beaver complex (ΔDOC = downstream − upstream)
   negligible, positive (DOC source) or negative (DOC sink)? → Figure 2
2. **Drivers** – Through which causal paths (hydrology: dams, dam height, channel gradient, water
   residence time; biology: phytoplankton, macrophytes, litter cover; solar radiation; incoming DOC)
   does beaver engineering change DOC? Structural causal model (SCM) fitted as linked Bayesian
   regressions → Figure 3
3. **Magnitude** – How large is the local change relative to the seasonal DOC contrast the stream
   already inherits from upstream, and does it reinforce, dampen or reverse that contrast?
   → Figure 4

## Repository structure

Everything that draws a figure or computes a number of the manuscript is code in `notebook/`;
every figure is written by the chunk that draws it **directly into `overleaf/`**, where the LaTeX
files read it. There is no copying step.

```
data/
  processed/            analysis data (described below)
  elevation/            Swiss elevation raster (map background of Figure 1)
images/                 hand-made illustrations that notebooks read (DAGs, schematic, photos)
notebook/
  manuscript/           Quarto documents: they draw the figures and compute the reported numbers
    analysis_plan.qmd       Figures 1–4 and SI Figures S7, S8, S11–S14 (read this first)
    appendix_s1_figures.qmd SI Figures S1, S3–S6 (main structural causal model)
    connectivity.qmd        SI Figure S10 (alternative model with lateral connectivity)
    concentrations_fit.qmd  per-stream concentration model that Figure 4 builds on
  fitting/              R scripts that fit the heavy brms models and cache their results
  references.bib, ecology-letters.csl   bibliography and citation style for the notebooks
results/                cached model summaries (*.rds); fitted models go to results/models_fit/ (not tracked)
overleaf/               LaTeX sources: main.tex, supplementary/supporting_information.tex and the
                        figures they include (figures/media/, supplementary/figures/media/)
```

## Data

All files are in `data/processed/`. DOC is in mg C L⁻¹.

### `m6.csv` — analysis dataset (360 rows)

The single analysis table, used by every notebook. One row per sampling (stream × season): 179
streams, upstream and downstream of the beaver complex (rows are paired by `site` and
`beaver_territory`). 158 streams have both seasons complete and a floodplain area; these are the
streams of the seasonal and magnitude analyses (those with `A_floodplain_km2`), while the
structural causal model and the direction model use all rows. 16 summer rows also carry the primary-producer
measurements; all other producer values are `NA` and are imputed in the models.

| Column(s) | Meaning (unit) |
|---|---|
| `site`, `beaver_territory`, `season` (`winter`/`summer`), `date`, `time` | Stream identifier, beaver territory, season and sampling time |
| `x_upstream`, `y_upstream`, `x_downstream`, `y_downstream` | Longitude/latitude of the upstream and downstream sampling points (WGS84) |
| `DOC_input`, `DOC_out`, `delta_DOC` | Upstream DOC, downstream DOC, and ΔDOC = `DOC_out` − `DOC_input` (mg L⁻¹) |
| `n_dams`, `dam_height`, `dam_persistence` | Dams upstream of the downstream point (count), maximum dam height (m), years of beaver occupation |
| `slope` | Channel gradient (%) from the DEM and the river line |
| `water_volume` | Pond water volume (m³) from imagery + DEM |
| `discharge` | Stream discharge at sampling from the PREVAH hydrological model (**L s⁻¹**) |
| `water_res_time` | Water residence time = `water_volume` / (3.6 × `discharge`) (**hours**) |
| `solar` | Solar radiation of the reach (Wh m⁻²) |
| `catchment_area_km` | Catchment area of the stream (km²), from the raw floodplain/catchment extraction (`catchArea`, m², divided by 10⁶). Not used by any model |
| `area_m6_revier` | Ponded (beaver-engineered) area in m². `area_m6_revier` / 10⁶ is `A_beaver_km2` |
| `plankton_abun` | Phytoplankton abundance (count of particles ≈10 µm–1 cm, from 50 L of pond water concentrated to 2 L through a 10 µm mesh and counted with a dark-field imaging microscope); 16 summer ponds, `NA` elsewhere |
| `macrophy_abun` | Macrophyte abundance (count of individual plants in the beaver pond); 16 summer ponds, `NA` elsewhere |
| `Cover_litter` | Soil litter cover (all dead plant material, including twigs < 7 cm circumference), estimated in a 1 × 5 m plot 0.5 m from the pond edge (% cover, 0–100 %); 16 summer ponds, `NA` elsewhere |
| `A_floodplain_km2` | Contributing floodplain area of the stream (km²): the area over which the inherited seasonal DOC contrast `S` is generated. Raw extraction (`fpArea`, m², divided by 10⁶) from the data provider; the same value in both seasons; `NA` for the 16 producer-measurement ponds (32 rows), which are not in the extraction. Rows with a floodplain area and DOC in both seasons form the 158 streams of the seasonal and magnitude analyses |
| `A_beaver_km2` | Ponded (beaver-engineered) area (km²) = `area_m6_revier` / 10⁶. The footprint contrast `F` = `A_floodplain_km2` / `A_beaver_km2` is computed in `analysis_plan.qmd` |
| `*_z` | The same variable standardised (mean 0, SD 1 over the analysed rows): `DOC_input_z`, `dam_height_z`, `slope_z`, `discharge_z`, `water_volume_z`, `water_res_time_z`, `solar_z`, `dam_persistence_z`, `n_dams_z`, `plankton_abun_z`, `cover_litter_z`, `macrophy_abun_z`. The models use the `_z` columns |

> **Units note.** Discharge is in L s⁻¹ and volume in m³, so the residence time is
> volume / (3.6 × discharge) in hours.

### Other files

- `later_connectivity_berger.csv` — `site`, `lateral_connectivity` (stream–wetland connectivity class: no wetland / wetland not connected / wetland connected; classification by K. Berger). Used by `connectivity.qmd`.
- `producer_sites_coordinates.csv` — `site`, `latitude`, `longitude` of the 16 producer-measurement sites (Figure 1a).

## How to reproduce the figures

**To reproduce a figure, render the notebook named in the tables below** (or run its chunk
interactively). The notebook writes the PNG to `overleaf/`, where `main.tex` or
`supplementary/supporting_information.tex` picks it up.

### Software

R ≥ 4.3, [Quarto](https://quarto.org) and [CmdStan](https://mc-stan.org/cmdstanr/)
(`cmdstanr::install_cmdstan()`), then:

```r
install.packages(c("brms","posterior","loo","tidyverse","tidybayes","ggdist","bayesplot",
  "dagitty","sensemakr","here","patchwork","ggpubr","ggh4x","ggpmisc","scales","knitr","DT",
  "RColorBrewer","NatParksPalettes","ggokabeito","modelr","terra","tidyterra","png","chron","readxl"))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
```

### Step 1 — fit the models (about an hour on 10 cores)

Fitted models are large and not tracked by git, so a fresh clone fits them once. From the
repository root:

```bash
Rscript notebook/fitting/me_refits.R                  # main SCM: delta_SCM_me, downstream_DOC_me
Rscript notebook/fitting/me_predictions.R             # predictions and leave-one-out
Rscript notebook/fitting/me_holdout.R                 # 80/20 hold-out by stream
Rscript notebook/fitting/dag_sensitivity_me.R         # alternative DAGs, prior sensitivity
Rscript notebook/fitting/downstream_DOC_me_fixed_fit.R  # inheritance slope fixed at 1
Rscript notebook/fitting/scm_downstream_fit.R         # plain SCM, input of the next script
Rscript notebook/fitting/pathway_decomposition.R      # pathway decomposition
quarto render notebook/manuscript/concentrations_fit.qmd   # per-stream concentrations (Figure 4)
```

Fits use `seed = 7` and 4 chains and are cached in `results/models_fit/` and `results/*.rds`;
rerunning a document reuses the cache.

### Step 2 — render the notebooks

```bash
quarto render notebook/manuscript/analysis_plan.qmd
quarto render notebook/manuscript/appendix_s1_figures.qmd
quarto render notebook/manuscript/connectivity.qmd
```

Models fitted inside the notebooks (the Figure 2 model, the intensity regressions, the negative
controls, the connectivity model) are cached the same way.

### Main-text figures

| Figure | File in `overleaf/figures/media/` | Notebook (section) | Needs |
|---|---|---|---|
| 1. Map of sites and schematic of the derived variables | `figure1_v2.png` | `analysis_plan.qmd` (11, Methods figure) | `m6.csv`, `producer_sites_coordinates.csv`, `data/elevation/`, `images/schematic_S_R_F.png` |
| 2. Direction of ΔDOC (posterior predictive, by season) | `Figure_2a.png` | `analysis_plan.qmd` (3, Q1) | `m6.csv` (the model is fitted in the chunk) |
| 3. Structural causal model of downstream DOC | `figure4_v1.png` | `analysis_plan.qmd` (4, Q2) | `results/me_refits.rds` (`me_refits.R`) |
| 4. Relative response *R<sub>s</sub>* and seasonal gain *G* | `figure3_v2.png` | `analysis_plan.qmd` (5, Q3) | `results/concentrations_thin.rds` (`concentrations_fit.qmd`) |

### Supporting-information figures and tables

Numbering is the order in `supporting_information.tex`; the file names keep an earlier numbering.
Files are in `overleaf/supplementary/figures/media/`.

| SI figure | File | Notebook (section) | Needs |
|---|---|---|---|
| S1 Stream-level ΔDOC under the main SCM | `figS7_site_level.png` | `appendix_s1_figures.qmd` | `me_refits.R` |
| S2 Hypothesised DAG (Table S1) | `figS1_dag_hypotheses.png` | hand-drawn, no code | — |
| S3 Priors | `figS2_priors.png` | `appendix_s1_figures.qmd` | — |
| S4 Posterior parameter distributions | `figS5_posteriors.png` | `appendix_s1_figures.qmd` | `me_refits.R` |
| S5 Residuals | `figS3_residuals.png` | `appendix_s1_figures.qmd` | `me_refits.R` |
| S6 Posterior predictive check | `figS4_ppc.png` | `appendix_s1_figures.qmd` | `me_refits.R` |
| S7 Implied conditional independencies | `figS10_dag_tests.png` | `analysis_plan.qmd` (9) | `m6.csv` |
| S8 Alternative DAGs, LOO and mediator effects | `figS11_dag_compare.png` | `analysis_plan.qmd` (9) | `dag_sensitivity_me.R` |
| S9 Extended DAG with lateral connectivity | `figS8_dag_connectivity_alt.png` | hand-drawn, no code | — |
| S10 Connectivity effects | `figS9_connectivity.png` | `connectivity.qmd` | `m6.csv`, `later_connectivity_berger.csv` |
| S11 Unmeasured-confounding robustness | `figS12_sensemakr.png` | `analysis_plan.qmd` (9) | `m6.csv` |
| S12 Negative controls | `figS13_negative_controls.png` | `analysis_plan.qmd` (9) | `m6.csv` |
| S13 Prior sensitivity | `figS14_prior_sens.png` | `analysis_plan.qmd` (9) | `dag_sensitivity_me.R` |
| S14 Spatial autocorrelation of residuals | `figS15_spatial_resid.png` | `analysis_plan.qmd` (9) | `me_predictions.R` |

Table S1 (hypotheses) is written by hand in the `.tex`; Table S2 (causal knowledge analysis) is built
in `analysis_plan.qmd` (`tbl-cka`). The numbers quoted in the text come from the same notebooks
(hold-out validation, intensity of beaver engineering, pathway decomposition: `analysis_plan.qmd`,
sections 6, 8 and 10).

### Building the PDFs

`overleaf/main.tex` and `overleaf/supplementary/supporting_information.tex` compile with any
standard LaTeX distribution (or in Overleaf, importing this repository); built PDFs are not tracked.

## License

GPL-3.0 (see `LICENSE`).
