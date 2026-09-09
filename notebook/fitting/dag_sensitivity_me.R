# Structural and prior sensitivity for the extended downstream-DOC SCM (downstream_DOC_me),
# which is now the manuscript's Q2 model: latent (measurement-error-corrected) upstream DOC, and
# the two causal links upstream DOC -> macrophytes, -> phytoplankton.
#
# DAG 0 is downstream_DOC_me itself (loaded, not refit). Three alternatives:
#   DAG 1 -- WITHOUT the two upstream -> producer links (the pre-extension structure): are the
#            mediator effects on downstream DOC an artefact of adding those links?
#   DAG 2 -- solar radiation acts directly on downstream DOC (photodegradation), not only through
#            the producers.
#   DAG 3 -- season is a common cause of the producers and of downstream DOC, left unblocked in
#            DAG 0.
# Plus prior sensitivity: every informative prior's sd halved and doubled, on DAG 0's structure.
# Sourced by analysis_plan.qmd, section 9. Cached in results/dag_sensitivity_me.rds.
suppressMessages({library(brms); library(dplyr); library(loo); library(here)})
m6 <- read.csv(here("data/processed/m6.csv")) |> mutate(across(c(season, site), as.factor))

med_bf <- function(up) list(
  bf(solar_z | mi() ~ 0 + Intercept, family = gaussian()),
  bf(as.formula(sprintf("plankton_abun_z | mi() ~ 0 + Intercept + mi(solar_z) + mi(water_res_time_z) + mi(%s)", up)), family = gaussian()),
  bf(as.formula(sprintf("macrophy_abun_z | mi() ~ 0 + Intercept + mi(solar_z) + mi(water_res_time_z) + mi(%s)", up)), family = gaussian()),
  bf(cover_litter_z | mi() ~ 0 + Intercept + mi(solar_z), family = gaussian()),
  bf(water_res_time_z | mi() ~ 0 + n_dams_z + mi(dam_height_z) + slope_z, family = gaussian()),
  bf(dam_height_z | mi() ~ 0 + Intercept, family = gaussian()))
med_bf_noup <- list(   # DAG 1: no upstream -> producers links
  bf(solar_z | mi() ~ 0 + Intercept, family = gaussian()),
  bf(plankton_abun_z | mi() ~ 0 + Intercept + mi(solar_z) + mi(water_res_time_z), family = gaussian()),
  bf(macrophy_abun_z | mi() ~ 0 + Intercept + mi(solar_z) + mi(water_res_time_z), family = gaussian()),
  bf(cover_litter_z | mi() ~ 0 + Intercept + mi(solar_z), family = gaussian()),
  bf(water_res_time_z | mi() ~ 0 + n_dams_z + mi(dam_height_z) + slope_z, family = gaussian()),
  bf(dam_height_z | mi() ~ 0 + Intercept, family = gaussian()))
add_season <- function(lst) { lst[[2]] <- update(lst[[2]], . ~ . + season); lst[[3]] <- update(lst[[3]], . ~ . + season)
  lst[[4]] <- update(lst[[4]], . ~ . + season); lst }
add_all <- function(x, lst) { for (b in lst) x <- x + b; x }

f_out <- bf(DOC_out | mi(0.2) ~ 0 + Intercept + mi(DOC_input) + mi(macrophy_abun_z) + mi(plankton_abun_z) +
              mi(cover_litter_z) + mi(water_res_time_z) + (0 + Intercept | site), family = student())
f_up  <- bf(DOC_input | mi(0.2) ~ 0 + Intercept, family = gaussian())

dags <- list(
  dag1 = add_all(f_out + f_up, med_bf_noup),
  dag2 = add_all(update(f_out, . ~ . + mi(solar_z)) + f_up, med_bf("DOC_input")),
  dag3 = add_all(update(f_out, . ~ . + season) + f_up, add_season(med_bf("DOC_input"))))
for (nm in names(dags)) dags[[nm]] <- dags[[nm]] + set_rescor(FALSE)

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
  set_prior("normal(0, 0.5)", class = "sigma", resp = c("solarz","coverlitterz","damheightz","macrophyabunz","planktonabunz","waterrestimez")),
  set_prior("normal(3, 3)", coef = "Intercept", resp = "DOCinput"), set_prior("student_t(3, 0, 2.5)", class = "sigma", resp = "DOCinput"),
  set_prior("normal(0, 2.5)", coef = "Intercept", resp = "DOCout"),
  set_prior("normal(1, 0.3)", coef = "miDOC_input", resp = "DOCout"),
  set_prior("normal(0.2, 0.075)", coef = c("micover_litter_z","mimacrophy_abun_z","miplankton_abun_z"), resp = "DOCout"),
  set_prior("normal(0, 0.1)", coef = "miwater_res_time_z", resp = "DOCout"),
  set_prior("gamma(2, 0.1)", class = "nu", lb = 2, resp = "DOCout"), set_prior("student_t(3, 0, 2.5)", class = "sigma", resp = "DOCout"),
  set_prior("normal(0, 2)", class = "sd", coef = "Intercept", group = "site", resp = "DOCout"))
pr_up_links <- set_prior("normal(0, 0.1)", coef = "miDOC_input", resp = c("planktonabunz", "macrophyabunz"))  # the two DAG-0 links

pr <- list(
  dag1 = common,                                                                     # no upstream -> producers links: drop their prior
  dag2 = c(common, pr_up_links, set_prior("normal(0, 0.5)", coef = "misolar_z", resp = "DOCout")),  # new direct path
  dag3 = c(common, pr_up_links, set_prior("normal(0, 0.5)", class = "b", resp = c("DOCout","planktonabunz","macrophyabunz","coverlitterz"))))  # season dummies

fit_one <- function(form, prior, tag) brm(form, data = m6, prior = prior, iter = 4000, warmup = 1500, chains = 4,
  cores = 10, threads = threading(2), control = list(max_treedepth = 12, adapt_delta = 0.9), backend = "cmdstanr", seed = 7,
  file = here("results", "models_fit", paste0("scm_me_", tag)), refresh = 0)

fits <- list(dag0 = readRDS(here("results", "models_fit", "downstream_DOC_me.rds")))
for (nm in names(dags)) { cat("fitting", nm, "\n"); fits[[nm]] <- fit_one(dags[[nm]], pr[[nm]], nm) }

# ---- prior sensitivity: halve and double the sd of every informative prior, on DAG 0's structure
pr_dag0 <- c(common, pr_up_links)
pr_half <- pr_dag0; pr_double <- pr_dag0
for (i in seq_len(nrow(pr_dag0))) {
  m <- regmatches(pr_dag0$prior[i], regexec("^normal\\(([-0-9.]+), ([0-9.]+)\\)$", pr_dag0$prior[i]))[[1]]
  if (length(m) == 3) { mu <- m[2]; s0 <- as.numeric(m[3])
    pr_half$prior[i] <- sprintf("normal(%s, %g)", mu, s0 / 2); pr_double$prior[i] <- sprintf("normal(%s, %g)", mu, s0 * 2) } }
form_dag0 <- add_all(f_out + f_up, med_bf("DOC_input")) + set_rescor(FALSE)
cat("fitting prior_half\n");   fits$prior_half   <- fit_one(form_dag0, pr_half, "prior_half")
cat("fitting prior_double\n"); fits$prior_double <- fit_one(form_dag0, pr_double, "prior_double")

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
saveRDS(list(models = summ, loo_compare = cmp), here("results", "dag_sensitivity_me.rds"))
print(cmp); cat("done\n")
