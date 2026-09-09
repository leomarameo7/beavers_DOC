# Both SCM parameterizations refitted with (i) measurement error on upstream
# DOC and (ii) the two new links upstream DOC -> phytoplankton, -> macrophytes.
# delta_SCM_me: dDOC outcome (original standardized upstream DOC, error 0.2/sd on the z scale)
# downstream_DOC_me: downstream-DOC outcome (upstream DOC in mg/L, error 0.2)
suppressMessages({library(brms); library(dplyr); library(here); library(posterior)})
m6 <- read.csv(here("data/processed/m6.csv")) |> mutate(across(c(season, site), as.factor))
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
  set_prior("normal(0, 0.1)", coef = "miDOC_input_z", resp = c("planktonabunz", "macrophyabunz")),   # the new links
  set_prior("normal(0, 2.5)", coef = "Intercept", resp = "deltaDOC"),
  set_prior("normal(0, 0.1)", coef = "miDOC_input_z", resp = "deltaDOC"),
  set_prior("normal(0.2, 0.075)", coef = c("micover_litter_z","mimacrophy_abun_z","miplankton_abun_z"), resp = "deltaDOC"),
  set_prior("normal(0, 0.1)", coef = "miwater_res_time_z", resp = "deltaDOC"),
  set_prior("gamma(2, 0.1)", class = "nu", lb = 2, resp = "deltaDOC"), set_prior("student_t(3, 0, 2.5)", class = "sigma", resp = "deltaDOC"),
  set_prior("normal(0, 2)", class = "sd", coef = "Intercept", group = "site", resp = "deltaDOC"))
pr_down <- c(common,
  set_prior("normal(3, 3)", coef = "Intercept", resp = "DOCinput"), set_prior("student_t(3, 0, 2.5)", class = "sigma", resp = "DOCinput"),
  set_prior("normal(0, 0.1)", coef = "miDOC_input", resp = c("planktonabunz", "macrophyabunz")),     # the new links
  set_prior("normal(0, 2.5)", coef = "Intercept", resp = "DOCout"),
  set_prior("normal(1, 0.3)", coef = "miDOC_input", resp = "DOCout"),
  set_prior("normal(0.2, 0.075)", coef = c("micover_litter_z","mimacrophy_abun_z","miplankton_abun_z"), resp = "DOCout"),
  set_prior("normal(0, 0.1)", coef = "miwater_res_time_z", resp = "DOCout"),
  set_prior("gamma(2, 0.1)", class = "nu", lb = 2, resp = "DOCout"), set_prior("student_t(3, 0, 2.5)", class = "sigma", resp = "DOCout"),
  set_prior("normal(0, 2)", class = "sd", coef = "Intercept", group = "site", resp = "DOCout"))
fit_one <- function(f, pr, tag) brm(f, data = m6, prior = pr, iter = 4000, warmup = 1500, chains = 4, cores = 10, threads = threading(2),
  control = list(max_treedepth = 12, adapt_delta = 0.9), backend = "cmdstanr", seed = 7, file = here("results","models_fit", tag), refresh = 0)
cat("fitting delta_SCM_me\n"); fd <- fit_one(f_delta, pr_delta, "delta_SCM_me")
cat("fitting downstream_DOC_me\n"); fo <- fit_one(f_down, pr_down, "downstream_DOC_me")
summ <- function(fit) { set.seed(7); dr <- as_draws_df(fit) |> subset_draws(draw = sort(sample(ndraws(fit), 2000)))
  nm <- grep("^bsp_|^b_", names(dr), value = TRUE); nm <- nm[!grepl("Intercept", nm)]
  list(tab = data.frame(par = nm, mean = sapply(nm, function(v) mean(dr[[v]])), lo = sapply(nm, function(v) quantile(dr[[v]], .025)),
                        hi = sapply(nm, function(v) quantile(dr[[v]], .975)), p_pos = sapply(nm, function(v) mean(dr[[v]] > 0))),
       draws = as.matrix(dr[, nm]), rhat = max(rhat(fit), na.rm = TRUE), ebfmi = NA) }
saveRDS(list(delta = summ(fd), down = summ(fo), s_up = s_up, se_z = se_z), here("results", "me_refits.rds"))
print(summ(fd)$tab[, -1] |> round(3)); print(summ(fo)$tab[, -1] |> round(3)); cat("done\n")
