# Robustness check on the extended Q2_downstream model (downstream_DOC_me): does fixing the
# upstream-DOC inheritance slope at 1 instead of estimating it freely change anything else?
# Same formula and priors as downstream_DOC_me (notebook/fitting/me_refits.R), with one change:
# the prior on the inheritance slope is tightened from Normal(1, 0.3) to Normal(1, 0.001), which
# pins it at 1 for practical purposes while keeping the mi() imputation of upstream DOC (an
# offset() term cannot carry missing/measurement-error data, so this is the standard way to fix
# an mi() coefficient in brms). Sourced by analysis_plan.qmd, section 4.
suppressMessages({library(brms); library(dplyr); library(here)})
m6 <- read.csv(here("data/processed/m6.csv")) |> mutate(across(c(season, site), as.factor))
s_up <- sd(m6$DOC_input, na.rm = TRUE)

med_bf <- function(up) list(
  bf(solar_z | mi() ~ 0 + Intercept, family = gaussian()),
  bf(as.formula(sprintf("plankton_abun_z | mi() ~ 0 + Intercept + mi(solar_z) + mi(water_res_time_z) + mi(%s)", up)), family = gaussian()),
  bf(as.formula(sprintf("macrophy_abun_z | mi() ~ 0 + Intercept + mi(solar_z) + mi(water_res_time_z) + mi(%s)", up)), family = gaussian()),
  bf(cover_litter_z | mi() ~ 0 + Intercept + mi(solar_z), family = gaussian()),
  bf(water_res_time_z | mi() ~ 0 + n_dams_z + mi(dam_height_z) + slope_z, family = gaussian()),
  bf(dam_height_z | mi() ~ 0 + Intercept, family = gaussian()))
add_all <- function(x, lst) { for (b in lst) x <- x + b; x }
f_down <- add_all(bf(DOC_out | mi(0.2) ~ 0 + Intercept + mi(DOC_input) + mi(macrophy_abun_z) + mi(plankton_abun_z) +
                       mi(cover_litter_z) + mi(water_res_time_z) + (0 + Intercept | site), family = student()) +
                  bf(DOC_input | mi(0.2) ~ 0 + Intercept, family = gaussian()), med_bf("DOC_input")) + set_rescor(FALSE)

pr_down_fixed <- c(
  set_prior("normal(0, 1)", coef = "Intercept", resp = c("coverlitterz","damheightz","macrophyabunz","planktonabunz","solarz")),
  set_prior("normal(-0.1, 0.075)", coef = "misolar_z", resp = "coverlitterz"),
  set_prior("normal(0.1, 0.075)",  coef = "misolar_z", resp = "macrophyabunz"),
  set_prior("normal(0, 0.075)",    coef = "miwater_res_time_z", resp = "macrophyabunz"),
  set_prior("normal(0.1, 0.075)",  coef = "misolar_z", resp = "planktonabunz"),
  set_prior("normal(0.1, 0.075)",  coef = "miwater_res_time_z", resp = "planktonabunz"),
  set_prior("normal(0.15, 0.075)", coef = "midam_height_z", resp = "waterrestimez"),
  set_prior("normal(0.15, 0.075)", coef = "n_dams_z", resp = "waterrestimez"),
  set_prior("normal(0, 0.075)",    coef = "slope_z", resp = "waterrestimez"),
  set_prior("normal(0, 0.5)", class = "sigma", resp = c("solarz","coverlitterz","damheightz","macrophyabunz","planktonabunz","waterrestimez")),
  set_prior("normal(3, 3)", coef = "Intercept", resp = "DOCinput"), set_prior("student_t(3, 0, 2.5)", class = "sigma", resp = "DOCinput"),
  set_prior("normal(0, 0.1)", coef = "miDOC_input", resp = c("planktonabunz", "macrophyabunz")),
  set_prior("normal(0, 2.5)", coef = "Intercept", resp = "DOCout"),
  set_prior("normal(1, 0.001)", coef = "miDOC_input", resp = "DOCout"),   # inheritance slope FIXED at 1
  set_prior("normal(0.2, 0.075)", coef = c("micover_litter_z","mimacrophy_abun_z","miplankton_abun_z"), resp = "DOCout"),
  set_prior("normal(0, 0.1)", coef = "miwater_res_time_z", resp = "DOCout"),
  set_prior("gamma(2, 0.1)", class = "nu", lb = 2, resp = "DOCout"), set_prior("student_t(3, 0, 2.5)", class = "sigma", resp = "DOCout"),
  set_prior("normal(0, 2)", class = "sd", coef = "Intercept", group = "site", resp = "DOCout"))

fit <- brm(f_down, data = m6, prior = pr_down_fixed, iter = 4000, warmup = 1500, chains = 4, cores = 10,
           threads = threading(2), control = list(max_treedepth = 12, adapt_delta = 0.9), backend = "cmdstanr",
           seed = 7, file = here("results", "models_fit", "downstream_DOC_me_fixed"), refresh = 0)
cat("DONE  max Rhat:", round(max(rhat(fit), na.rm = TRUE), 4), "\n")
d <- as_draws_df(fit)
for (v in c("bsp_DOCout_miDOC_input","bsp_DOCout_mimacrophy_abun_z","bsp_DOCout_miplankton_abun_z",
            "bsp_DOCout_micover_litter_z","bsp_DOCout_miwater_res_time_z","b_DOCout_Intercept"))
  cat(sprintf("%-34s %7.3f [%7.3f, %7.3f]\n", sub("bsp_DOCout_mi|b_DOCout_","",v), mean(d[[v]]), quantile(d[[v]],.025), quantile(d[[v]],.975)))
