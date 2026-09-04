# Posterior DOC concentrations per stream, used to derive S, B, G and R_s.
#
# What this script does, step by step:
#   1. Pairs the winter and summer sampling of each stream (158 streams have both).
#   2. Fits one Bayesian hierarchical model per measured concentration
#      (upstream winter, downstream winter, upstream summer, downstream summer),
#      on the log scale (concentrations are positive and right-skewed),
#      with the 0.2 mg/L instrument error attached to every observation
#      (transferred to the log scale by the delta method: se_log = 0.2 / value).
#   3. Extracts, for every posterior draw, each stream's estimated true
#      concentration, and computes the derived quantities from those draws —
#      so no ratio or difference is ever formed from raw values, and the
#      measurement uncertainty is carried into S, B, G and R_s.
# Output: results/concentrations.rds (full draws; a thinned copy
# results/concentrations_thin.rds with 4000 draws is tracked in git).

suppressMessages({library(brms); library(tidyverse); library(here)})

# ---- 1. One row per stream, with the four measured concentrations ----------
m <- read.csv(here("data", "processed", "m12.csv"))
w <- m |> filter(season == "winter") |>
  select(site, beaver_territory, up_w = DOC_input, dn_w = DOC_out,
         Afp = A_floodplain_km2, Abv = A_beaver_km2)          # winter sampling
s <- m |> filter(season == "summer") |>
  select(site, beaver_territory, up_s = DOC_input, dn_s = DOC_out)  # summer sampling
# keep only streams sampled in BOTH seasons with all four concentrations: 158 of 176
d <- inner_join(w, s, by = c("site", "beaver_territory")) |>
  filter(complete.cases(up_w, dn_w, up_s, dn_s)) |>
  mutate(reach = factor(row_number()))   # stream id used as the grouping factor

# ---- 2. One hierarchical measurement-error model per concentration ---------
# weakly informative priors: log-concentration centred near 1 (about 2.7 mg/L)
pri <- c(prior(normal(1, 1.5), class = "Intercept"),
         prior(student_t(3, 0, 1), class = "sd"))
fitc <- function(v) {
  dd <- d |> mutate(
    y      = log(pmax(.data[[v]], 0.05)),        # log concentration (floored at 0.05)
    se_log = 0.2 / pmax(.data[[v]], 0.25))       # instrument error on the log scale
  # y | se(...): each observation carries its known measurement error;
  # (1 | reach): one true-concentration deviation per stream
  brm(bf(y | se(se_log, sigma = FALSE) ~ 1 + (1 | reach)),
      data = dd, prior = pri, chains = 4, cores = 4, iter = 6000, warmup = 2000,
      seed = 7, backend = "cmdstanr",
      control = list(adapt_delta = 0.99, max_treedepth = 12),
      file = here("results", "models_fit", paste0("conc_", v)), refresh = 0)
}

# ---- 3. Per-draw true concentrations and the derived quantities ------------
# for each draw: intercept + stream deviation, back-transformed to mg/L
conc <- function(f) { x <- as_draws_df(f)
  exp(as.matrix(x[, grep("^r_reach\\[", names(x))]) + x$b_Intercept) }
L <- setNames(lapply(c("up_w", "dn_w", "up_s", "dn_s"), function(v) conc(fitc(v))),
              c("up_w", "dn_w", "up_s", "dn_s"))
S  <- L$up_s - L$up_w   # inherited seasonal contrast (upstream, summer - winter)
dw <- L$dn_w - L$up_w   # local response in winter
ds <- L$dn_s - L$up_s   # local response in summer
saveRDS(list(S = S, dw = dw, ds = ds,
             R_s = ds / abs(S), R_w = dw / abs(S),   # relative response by season
             up_s = L$up_s, up_w = L$up_w, d = d),
        here("results", "concentrations.rds"))
