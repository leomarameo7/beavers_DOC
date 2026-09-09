# Hold-out test of the two extended SCM parameterizations (Q2_delta =
# delta_SCM_me, Q2_downstream = downstream_DOC_me): refit both on 80% of the
# streams (same split as the original scm_holdout.R, seed 42) and predict at
# the other 20%, which neither model has seen -- no site effect, mediators
# mostly unmeasured there too. Complements the leave-one-out test already in
# results/me_predictions.rds. Sourced by scm_outcome_choice.qmd.
suppressMessages({library(brms); library(dplyr); library(here); library(posterior)})
m6 <- read.csv(here("data/processed/m6.csv")) |> mutate(across(c(season, site), as.factor))
set.seed(42)
sites <- unique(m6$site)
test_sites <- sample(sites, size = round(0.2 * length(sites)))
train <- m6 |> filter(!site %in% test_sites) |> mutate(site = droplevels(site))
test  <- m6 |> filter(site %in% test_sites, !is.na(DOC_out), !is.na(DOC_input))

s_up <- sd(m6$DOC_input, na.rm = TRUE); se_z <- round(0.2 / s_up, 4)

med_bf <- function(up) list(
  bf(solar_z | mi() ~ 0 + Intercept, family = gaussian()),
  bf(as.formula(sprintf("plankton_abun_z | mi() ~ 0 + Intercept + mi(solar_z) + mi(water_res_time_z) + mi(%s)", up)), family = gaussian()),
  bf(as.formula(sprintf("macrophy_abun_z | mi() ~ 0 + Intercept + mi(solar_z) + mi(water_res_time_z) + mi(%s)", up)), family = gaussian()),
  bf(cover_litter_z | mi() ~ 0 + Intercept + mi(solar_z), family = gaussian()),
  bf(water_res_time_z | mi() ~ 0 + n_dams_z + mi(dam_height_z) + slope_z, family = gaussian()),
  bf(dam_height_z | mi() ~ 0 + Intercept, family = gaussian()))
add_all <- function(x, lst) { for (b in lst) x <- x + b; x }

f_delta <- add_all(bf(delta_DOC | mi(0.283) ~ 0 + Intercept + mi(DOC_input_z) + mi(macrophy_abun_z) + mi(plankton_abun_z) +
                        mi(cover_litter_z) + mi(water_res_time_z) + (0 + Intercept | site), family = student()) +
                   bf(as.formula(sprintf("DOC_input_z | mi(%s) ~ 0 + Intercept", se_z)), family = gaussian()), med_bf("DOC_input_z")) + set_rescor(FALSE)
f_down  <- add_all(bf(DOC_out | mi(0.2) ~ 0 + Intercept + mi(DOC_input) + mi(macrophy_abun_z) + mi(plankton_abun_z) +
                        mi(cover_litter_z) + mi(water_res_time_z) + (0 + Intercept | site), family = student()) +
                   bf(DOC_input | mi(0.2) ~ 0 + Intercept, family = gaussian()), med_bf("DOC_input")) + set_rescor(FALSE)

common <- c(
  set_prior("normal(0, 1)", coef = "Intercept", resp = c("coverlitterz","damheightz","macrophyabunz","planktonabunz","solarz")),
  set_prior("normal(-0.1, 0.075)", coef = "misolar_z", resp = "coverlitterz"),
  set_prior("normal(0.1, 0.075)",  coef = "misolar_z", resp = "macrophyabunz"),
  set_prior("normal(0, 0.075)",    coef = "miwater_res_time_z", resp = "macrophyabunz"),
  set_prior("normal(0.1, 0.075)",  coef = "misolar_z", resp = "planktonabunz"),
  set_prior("normal(0.1, 0.075)",  coef = "miwater_res_time_z", resp = "planktonabunz"),
  set_prior("normal(0.15, 0.075)", coef = "midam_height_z", resp = "waterrestimez"),
  set_prior("normal(0.15, 0.075)", coef = "n_dams_z", resp = "waterrestimez"),
  set_prior("normal(0, 0.075)",    coef = "slope_z", resp = "waterrestimez"),
  set_prior("normal(0, 0.5)", class = "sigma", resp = c("solarz","coverlitterz","damheightz","macrophyabunz","planktonabunz","waterrestimez")))
pr_delta <- c(common,
  set_prior("normal(0, 1)", coef = "Intercept", resp = "DOCinputz"), set_prior("normal(0, 0.5)", class = "sigma", resp = "DOCinputz"),
  set_prior("normal(0, 0.1)", coef = "miDOC_input_z", resp = c("planktonabunz", "macrophyabunz")),
  set_prior("normal(0, 2.5)", coef = "Intercept", resp = "deltaDOC"),
  set_prior("normal(0, 0.1)", coef = "miDOC_input_z", resp = "deltaDOC"),
  set_prior("normal(0.2, 0.075)", coef = c("micover_litter_z","mimacrophy_abun_z","miplankton_abun_z"), resp = "deltaDOC"),
  set_prior("normal(0, 0.1)", coef = "miwater_res_time_z", resp = "deltaDOC"),
  set_prior("gamma(2, 0.1)", class = "nu", lb = 2, resp = "deltaDOC"), set_prior("student_t(3, 0, 2.5)", class = "sigma", resp = "deltaDOC"),
  set_prior("normal(0, 2)", class = "sd", coef = "Intercept", group = "site", resp = "deltaDOC"))
pr_down <- c(common,
  set_prior("normal(3, 3)", coef = "Intercept", resp = "DOCinput"), set_prior("student_t(3, 0, 2.5)", class = "sigma", resp = "DOCinput"),
  set_prior("normal(0, 0.1)", coef = "miDOC_input", resp = c("planktonabunz", "macrophyabunz")),
  set_prior("normal(0, 2.5)", coef = "Intercept", resp = "DOCout"),
  set_prior("normal(1, 0.3)", coef = "miDOC_input", resp = "DOCout"),
  set_prior("normal(0.2, 0.075)", coef = c("micover_litter_z","mimacrophy_abun_z","miplankton_abun_z"), resp = "DOCout"),
  set_prior("normal(0, 0.1)", coef = "miwater_res_time_z", resp = "DOCout"),
  set_prior("gamma(2, 0.1)", class = "nu", lb = 2, resp = "DOCout"), set_prior("student_t(3, 0, 2.5)", class = "sigma", resp = "DOCout"),
  set_prior("normal(0, 2)", class = "sd", coef = "Intercept", group = "site", resp = "DOCout"))

fit_one <- function(f, pr, tag) brm(f, data = train, prior = pr, iter = 4000, warmup = 1500, chains = 4, cores = 10, threads = threading(2),
  control = list(max_treedepth = 12, adapt_delta = 0.9), backend = "cmdstanr", seed = 7, file = here("results","models_fit", tag), refresh = 0)
cat("fitting delta_SCM_me on", nrow(train), "training rows (", length(unique(train$site)), "streams )\n")
fd <- fit_one(f_delta, pr_delta, "delta_SCM_me_train")
cat("fitting downstream_DOC_me on", nrow(train), "training rows\n")
fo <- fit_one(f_down, pr_down, "downstream_DOC_me_train")

# unmeasured-at-a-new-stream covariates set to their standardized mean (0), as in scm_holdout.R
zero_na <- function(d) d |> mutate(across(c(macrophy_abun_z, plankton_abun_z, cover_litter_z, water_res_time_z, solar_z, dam_height_z), ~ ifelse(is.na(.x), 0, .x)))
nd_test  <- zero_na(test)
nd_train <- zero_na(train |> filter(!is.na(DOC_out), !is.na(DOC_input)))

metrics <- function(y, mu, q = NULL) c(R2 = 1 - sum((y - mu)^2) / sum((y - mean(y))^2), RMSE = sqrt(mean((y - mu)^2)), RMSE_null = sd(y),
                                        cover90 = if (is.null(q)) NA else mean(y >= q[1, ] & y <= q[2, ]))

# Q2_delta: outcome dDOC, out-of-sample (re_formula = NA: unseen streams) and in-sample
ep_d_test  <- posterior_epred(fd,  resp = "deltaDOC", newdata = nd_test,  re_formula = NA)
pp_d_test  <- posterior_predict(fd, resp = "deltaDOC", newdata = nd_test,  re_formula = NA)
ep_d_train <- posterior_epred(fd,  resp = "deltaDOC", newdata = nd_train)
pp_d_train <- posterior_predict(fd, resp = "deltaDOC", newdata = nd_train)
mu_d_test  <- colMeans(ep_d_test);  q_d_test  <- apply(pp_d_test,  2, quantile, c(.05, .95))
mu_d_train <- colMeans(ep_d_train); q_d_train <- apply(pp_d_train, 2, quantile, c(.05, .95))
D_out <- metrics(nd_test$delta_DOC,  mu_d_test,  q_d_test)
D_in  <- metrics(nd_train$delta_DOC, mu_d_train, q_d_train)

# Q2_downstream: outcome DOC_out, and derived dDOC = predicted DOC_out - observed DOC_input
ep_o_test  <- posterior_epred(fo,  resp = "DOCout", newdata = nd_test,  re_formula = NA)
pp_o_test  <- posterior_predict(fo, resp = "DOCout", newdata = nd_test,  re_formula = NA)
ep_o_train <- posterior_epred(fo,  resp = "DOCout", newdata = nd_train)
pp_o_train <- posterior_predict(fo, resp = "DOCout", newdata = nd_train)
mu_o_test  <- colMeans(ep_o_test);  q_o_test  <- apply(pp_o_test,  2, quantile, c(.05, .95))
mu_o_train <- colMeans(ep_o_train); q_o_train <- apply(pp_o_train, 2, quantile, c(.05, .95))
O_out <- metrics(nd_test$DOC_out,  mu_o_test,  q_o_test)
O_in  <- metrics(nd_train$DOC_out, mu_o_train, q_o_train)

dd_test_pred  <- sweep(ep_o_test,  2, nd_test$DOC_input);  dd_test_obs  <- nd_test$DOC_out  - nd_test$DOC_input
dd_train_pred <- sweep(ep_o_train, 2, nd_train$DOC_input); dd_train_obs <- nd_train$DOC_out - nd_train$DOC_input
Od_out <- metrics(dd_test_obs,  colMeans(dd_test_pred))
Od_in  <- metrics(dd_train_obs, colMeans(dd_train_pred))
# predictive interval for the derived dDOC at held-out streams: unlike LOO, this is a genuine
# posterior predictive draw for unseen data, so no importance-sampling reweighting is needed.
dd_q_test <- apply(sweep(pp_o_test, 2, nd_test$DOC_input), 2, quantile, c(.05, .95))

saveRDS(list(test_sites = test_sites, n_test_streams = length(test_sites), n_test_rows = nrow(nd_test), n_train_rows = nrow(nd_train),
             rhat = c(delta = max(rhat(fd), na.rm = TRUE), down = max(rhat(fo), na.rm = TRUE)),
             D_out = D_out, D_in = D_in, O_out = O_out, O_in = O_in, Od_out = Od_out, Od_in = Od_in,
             obs_d_test = nd_test$delta_DOC, mu_d_test = mu_d_test, q_d_test = q_d_test,
             obs_o_test = nd_test$DOC_out, mu_o_test = mu_o_test, q_o_test = q_o_test,
             up_test = nd_test$DOC_input, dd_test_pred_med = apply(dd_test_pred, 2, median), dd_q_test = dd_q_test),
        here("results", "me_holdout.rds"))
cat("done. max Rhat:", round(max(max(rhat(fd), na.rm = TRUE), max(rhat(fo), na.rm = TRUE)), 3),
    " test streams:", length(test_sites), " test rows:", nrow(nd_test), "\n")
print(round(rbind(D_in, D_out, O_in, O_out, Od_in, Od_out), 3))
