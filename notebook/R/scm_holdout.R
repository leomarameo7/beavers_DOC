# Hold-out validation of the downstream-DOC SCM: fit on 80% of the streams,
# predict downstream DOC at the other 20%. Sourced by analysis_plan.qmd
# (chunk 'holdout-fit'); brms caches the fit in results/models_fit/.
suppressMessages({library(brms); library(dplyr); library(here)})
m6 <- read.csv(here("data/processed/m6.csv")) |> mutate(across(c(season, site), as.factor))
set.seed(42)
sites <- unique(m6$site)
test_sites <- sample(sites, size = round(0.2 * length(sites)))
train <- m6 |> filter(!site %in% test_sites)
test  <- m6 |> filter(site %in% test_sites, !is.na(DOC_out), !is.na(DOC_input))
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
  set_prior("normal(3, 3)", coef = "Intercept", resp = "DOCinput"),
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
  set_prior("normal(1, 0.3)",      coef = "miDOC_input", resp = "DOCout"),
  set_prior("normal(0.2, 0.075)",  coef = "micover_litter_z", resp = "DOCout"),
  set_prior("normal(0.2, 0.075)",  coef = "mimacrophy_abun_z", resp = "DOCout"),
  set_prior("normal(0.2, 0.075)",  coef = "miplankton_abun_z", resp = "DOCout"),
  set_prior("normal(0, 0.1)",      coef = "miwater_res_time_z", resp = "DOCout"),
  set_prior("gamma(2, 0.1)", class = "nu", lb = 2, resp = "DOCout"),
  set_prior("normal(0, 2)", class = "sd", coef = "Intercept", group = "site", resp = "DOCout"),
  set_prior("normal(0, 0.5)", class = "sigma", resp = c("solarz","coverlitterz","damheightz","macrophyabunz","planktonabunz","waterrestimez")))
fit_train <- brm(bform_down, data = train, prior = pr_down, iter = 4000, warmup = 1500, chains = 4, cores = 10,
                 threads = threading(2), control = list(max_treedepth = 12), backend = "cmdstanr", seed = 7,
                 file = here("results","models_fit","b1_downstream_train"), refresh = 0)
# prediction at unseen streams: mediators not measured there enter at their
# (standardized) mean of 0; site effects are unknown -> population-level prediction
nd <- test |> mutate(across(c(macrophy_abun_z, plankton_abun_z, cover_litter_z, water_res_time_z,
                             solar_z, dam_height_z), ~ ifelse(is.na(.x), 0, .x)))
ep <- posterior_epred(fit_train, resp = "DOCout", newdata = nd, re_formula = NA)
pp <- posterior_predict(fit_train, resp = "DOCout", newdata = nd, re_formula = NA)
# in-sample predictions for the training streams (site effects known)
tr <- train |> filter(!is.na(DOC_out), !is.na(DOC_input)) |> mutate(site = droplevels(site))
ep_tr <- posterior_epred(fit_train, resp = "DOCout", newdata = tr |> mutate(across(c(macrophy_abun_z, plankton_abun_z,
            cover_litter_z, water_res_time_z, solar_z, dam_height_z), ~ ifelse(is.na(.x), 0, .x))))
pp_tr <- posterior_predict(fit_train, resp = "DOCout", newdata = tr |> mutate(across(c(macrophy_abun_z, plankton_abun_z,
            cover_litter_z, water_res_time_z, solar_z, dam_height_z), ~ ifelse(is.na(.x), 0, .x))))
saveRDS(list(test = nd, epred = ep, pred = pp, test_sites = test_sites,
             train = tr, epred_train = ep_tr, pred_train = pp_tr,
             rhat = max(rhat(fit_train), na.rm = TRUE)), here("results", "holdout_pred.rds"))
cat("holdout done, max Rhat:", round(max(rhat(fit_train), na.rm = TRUE), 3), " test rows:", nrow(nd), "\n")
