suppressMessages({library(brms);library(dplyr);library(here)})
setwd("/Users/leonardocapitani/Documents/Git/beavers_DOC")
m6 <- read.csv("data/processed/m6.csv") |> mutate(across(c(season, site), as.factor))
# Robustness check requested by Josh (4 Sep 2026, 18:12): does fixing the
# upstream-DOC inheritance slope at 1 change any other path in the SCM?
# Same model and priors as b1_downstream (scm_downstream_fit.R), with one
# change: the prior on the inheritance slope is tightened from Normal(1, 0.3)
# to Normal(1, 0.001), which pins it at 1 for practical purposes while keeping
# the mi() imputation of the 20 missing upstream-DOC values (an offset() term
# cannot carry missing data, so this is the standard way to fix an mi()
# coefficient in brms).
bform_down <-
  bf(DOC_out | mi(0.2) ~ 0 + Intercept + mi(DOC_input) + mi(macrophy_abun_z) + mi(plankton_abun_z) +
       mi(cover_litter_z) + mi(water_res_time_z) + (0 + Intercept | site), family = student()) +
  bf(solar_z | mi() ~ 0 + Intercept, family = gaussian()) +
  bf(DOC_input | mi() ~ 0 + Intercept, family = gaussian()) +
  bf(plankton_abun_z | mi() ~ 0 + Intercept + mi(solar_z) + mi(water_res_time_z), family = gaussian()) +
  bf(macrophy_abun_z | mi() ~ 0 + Intercept + mi(solar_z) + mi(water_res_time_z), family = gaussian()) +
  bf(cover_litter_z | mi() ~ 0 + Intercept + mi(solar_z), family = gaussian()) +
  bf(water_res_time_z | mi() ~ 0 + n_dams_z + mi(dam_height_z) + slope_z, family = gaussian()) +
  bf(dam_height_z | mi() ~ 0 + Intercept, family = gaussian()) + set_rescor(FALSE)
pr_down <- c(
  set_prior("normal(0, 1)", coef = "Intercept", resp = c("coverlitterz","damheightz","macrophyabunz","planktonabunz","solarz")),
  set_prior("normal(3, 3)", coef = "Intercept", resp = "DOCinput"),          # raw mg/L scale
  set_prior("student_t(3, 0, 2.5)", class = "sigma", resp = "DOCinput"),
  set_prior("normal(-0.1, 0.075)", coef = "misolar_z", resp = "coverlitterz"),
  set_prior("normal(0.1, 0.075)",  coef = "misolar_z", resp = "macrophyabunz"),
  set_prior("normal(0, 0.075)",    coef = "miwater_res_time_z", resp = "macrophyabunz"),
  set_prior("normal(0.1, 0.075)",  coef = "misolar_z", resp = "planktonabunz"),
  set_prior("normal(0.1, 0.075)",  coef = "miwater_res_time_z", resp = "planktonabunz"),
  set_prior("normal(0.15, 0.075)", coef = "midam_height_z", resp = "waterrestimez"),
  set_prior("normal(0.15, 0.075)", coef = "n_dams_z", resp = "waterrestimez"),
  set_prior("normal(0, 0.075)",    coef = "slope_z", resp = "waterrestimez"),
  set_prior("normal(0, 2.5)",      coef = "Intercept", resp = "DOCout"),
  set_prior("normal(1, 0.001)",    coef = "miDOC_input", resp = "DOCout"),   # inheritance slope FIXED at 1
  set_prior("normal(0.2, 0.075)",  coef = "micover_litter_z", resp = "DOCout"),
  set_prior("normal(0.2, 0.075)",  coef = "mimacrophy_abun_z", resp = "DOCout"),
  set_prior("normal(0.2, 0.075)",  coef = "miplankton_abun_z", resp = "DOCout"),
  set_prior("normal(0, 0.1)",      coef = "miwater_res_time_z", resp = "DOCout"),
  set_prior("gamma(2, 0.1)", class = "nu", lb = 2, resp = "DOCout"),
  set_prior("normal(0, 2)", class = "sd", coef = "Intercept", group = "site", resp = "DOCout"),
  set_prior("normal(0, 0.5)", class = "sigma", resp = c("solarz","coverlitterz","damheightz","macrophyabunz","planktonabunz","waterrestimez")))
fit <- brm(bform_down, data = m6, prior = pr_down, iter = 10000, warmup = 3000, chains = 4, cores = 10,
           threads = threading(2), control = list(max_treedepth = 12), backend = "cmdstanr", seed = 7,
           file = here("results","models_fit","b1_downstream_fixed"), refresh = 0)
cat("DONE  max Rhat:", round(max(rhat(fit), na.rm = TRUE), 4), "\n")
d <- as_draws_df(fit)
for (v in c("bsp_DOCout_miDOC_input","bsp_DOCout_mimacrophy_abun_z","bsp_DOCout_miplankton_abun_z",
            "bsp_DOCout_micover_litter_z","bsp_DOCout_miwater_res_time_z","b_DOCout_Intercept"))
  cat(sprintf("%-34s %7.3f [%7.3f, %7.3f]\n", sub("bsp_DOCout_mi|b_DOCout_","",v), mean(d[[v]]), quantile(d[[v]],.025), quantile(d[[v]],.975)))
