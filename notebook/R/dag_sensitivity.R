# Structural sensitivity (alternative DAGs) and prior sensitivity for the
# downstream-DOC SCM. Sourced by analysis_plan.qmd; fits are cached by brms in
# results/models_fit/, summaries in results/dag_sensitivity.rds.
suppressMessages({library(brms); library(dplyr); library(loo); library(here)})
m6 <- read.csv(here("data/processed/m6.csv")) |> mutate(across(c(season, site), as.factor))

# ---- the current model (DAG 0), exactly as b1_downstream ---------------------
f_out  <- DOC_out | mi(0.2) ~ 0 + Intercept + mi(DOC_input) + mi(macrophy_abun_z) + mi(plankton_abun_z) +
            mi(cover_litter_z) + mi(water_res_time_z) + (0 + Intercept | site)
f_pla  <- plankton_abun_z | mi() ~ 0 + Intercept + mi(solar_z) + mi(water_res_time_z)
f_mac  <- macrophy_abun_z | mi() ~ 0 + Intercept + mi(solar_z) + mi(water_res_time_z)
f_lit  <- cover_litter_z  | mi() ~ 0 + Intercept + mi(solar_z)
build <- function(fo, fp, fm, fl) bf(fo, family = student()) + bf(fp, family = gaussian()) +
  bf(fm, family = gaussian()) + bf(fl, family = gaussian()) +
  bf(solar_z | mi() ~ 0 + Intercept, family = gaussian()) +
  bf(DOC_input | mi() ~ 0 + Intercept, family = gaussian()) +
  bf(water_res_time_z | mi() ~ 0 + n_dams_z + mi(dam_height_z) + slope_z, family = gaussian()) +
  bf(dam_height_z | mi() ~ 0 + Intercept, family = gaussian()) + set_rescor(FALSE)

pr_base <- c(
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
# class-level fallbacks for the paths added by the alternative DAGs (coef priors above take precedence)
pr_new <- c(set_prior("normal(0, 0.1)", class = "b", resp = c("planktonabunz","macrophyabunz","coverlitterz")),
            set_prior("normal(0, 0.5)", class = "b", resp = "DOCout"))

dags <- list(
  dag0 = build(f_out, f_pla, f_mac, f_lit),
  # DAG 1: upstream DOC also drives the primary producers (Leonardo's alternative)
  dag1 = build(f_out,
               update(f_pla, . ~ . + mi(DOC_input)), update(f_mac, . ~ . + mi(DOC_input)), f_lit),
  # DAG 2: solar radiation acts directly on downstream DOC as well as through the producers
  dag2 = build(update(f_out, . ~ . + mi(solar_z)), f_pla, f_mac, f_lit),
  # DAG 3: season is a common cause of the producers and of downstream DOC (unblocked confounder)
  dag3 = build(update(f_out, . ~ . + season), update(f_pla, . ~ . + season),
               update(f_mac, . ~ . + season), update(f_lit, . ~ . + season)))

fit_one <- function(form, prior, tag) brm(form, data = m6, prior = prior, iter = 4000, warmup = 1500, chains = 4,
  cores = 10, threads = threading(2), control = list(max_treedepth = 12), backend = "cmdstanr", seed = 7,
  file = here("results", "models_fit", paste0("scm_", tag)), refresh = 0)

fits <- list()
for (nm in names(dags)) { cat("fitting", nm, "\n"); fits[[nm]] <- fit_one(dags[[nm]], c(pr_base, pr_new), nm) }

# ---- prior sensitivity: halve and double the sd of every informative prior ----
pr_all <- c(pr_base, pr_new)
pr_half <- pr_all; pr_double <- pr_all
for (i in seq_len(nrow(pr_all))) {
  m <- regmatches(pr_all$prior[i], regexec("^normal\\(([-0-9.]+), ([0-9.]+)\\)$", pr_all$prior[i]))[[1]]
  if (length(m) == 3) { mu <- m[2]; s0 <- as.numeric(m[3])
    pr_half$prior[i] <- sprintf("normal(%s, %g)", mu, s0 / 2); pr_double$prior[i] <- sprintf("normal(%s, %g)", mu, s0 * 2) } }
cat("fitting prior_half\n");   fits$prior_half   <- fit_one(dags$dag0, pr_half, "prior_half")
cat("fitting prior_double\n"); fits$prior_double <- fit_one(dags$dag0, pr_double, "prior_double")

# ---- summaries: mediator effects, new paths, LOO on downstream DOC ------------
obs <- which(!is.na(m6$DOC_out))
keep <- c("bsp_DOCout_miDOC_input", "bsp_DOCout_mimacrophy_abun_z", "bsp_DOCout_miplankton_abun_z",
          "bsp_DOCout_micover_litter_z", "bsp_DOCout_miwater_res_time_z")
summ <- lapply(names(fits), function(nm) {
  f <- fits[[nm]]; d <- as_draws_df(f); nd <- ndraws(f)
  ll <- log_lik(f, resp = "DOCout")[, obs]
  lo <- loo(ll, r_eff = relative_eff(exp(ll), chain_id = rep(1:4, each = nd / 4)))
  extra <- setdiff(grep("^bsp_|^b_", names(d), value = TRUE), keep)
  extra <- extra[grepl("DOCout|planktonabunz|macrophyabunz|coverlitterz", extra) & !grepl("Intercept", extra)]
  list(name = nm, draws = as.matrix(d[sample(nd, 2000), keep]), extra = as.matrix(d[sample(nd, 2000), extra, drop = FALSE]),
       loo = list(estimates = lo$estimates, pareto_k = lo$diagnostics$pareto_k), loo_obj = lo,
       rhat = max(rhat(f), na.rm = TRUE), formula = format(f$formula))
})
names(summ) <- names(fits)
cmp <- loo_compare(lapply(summ[c("dag0","dag1","dag2","dag3")], `[[`, "loo_obj"))
summ <- lapply(summ, function(x) { x$loo_obj <- NULL; x })
saveRDS(list(models = summ, loo_compare = cmp), here("results", "dag_sensitivity.rds"))
print(cmp); cat("done\n")
