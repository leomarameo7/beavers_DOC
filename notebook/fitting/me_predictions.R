# Everything scm_outcome_choice.qmd needs from the two extended models
# (delta_SCM_me, downstream_DOC_me): comparison draws, predictions, LOO for each
# model's own outcome and for dDOC, derived dDOC. Cached in results/me_predictions.rds.
suppressMessages({library(brms); library(dplyr); library(loo); library(here); library(posterior)})
m6 <- read.csv(here("data/processed/m6.csv"))
fd <- readRDS(here("results/models_fit/delta_SCM_me.rds")); fo <- readRDS(here("results/models_fit/downstream_DOC_me.rds"))
stopifnot(nrow(fd$data) == nrow(m6), nrow(fo$data) == nrow(m6))
s_up <- sd(m6$DOC_input, na.rm = TRUE); mu_up <- mean(m6$DOC_input, na.rm = TRUE)
set.seed(7); nd <- ndraws(fo); keep <- sort(sample(nd, 2000))
pick <- function(fit, pat) { d <- as_draws_df(fit) |> subset_draws(draw = keep); as.matrix(d[, grep(pat, names(d), value = TRUE)]) }
cmpD <- pick(fd, "^(b_deltaDOC_Intercept|bsp_deltaDOC_|bsp_planktonabunz_miDOC|bsp_macrophyabunz_miDOC|sd_site__deltaDOC|sigma_deltaDOC|nu_deltaDOC)")
cmpO <- pick(fo, "^(b_DOCout_Intercept|bsp_DOCout_|bsp_planktonabunz_miDOC|bsp_macrophyabunz_miDOC|sd_site__DOCout|sigma_DOCout|nu_DOCout)")
chain_id <- rep(1:4, each = nd / 4)
loo_block <- function(fit, resp, y, obs) {
  ep <- posterior_epred(fit, resp = resp); pp <- posterior_predict(fit, resp = resp); ll <- log_lik(fit, resp = resp)
  r_eff <- relative_eff(exp(ll[, obs]), chain_id = chain_id); lo <- loo(ll[, obs], r_eff = r_eff, save_psis = TRUE)
  mu <- E_loo(ep[, obs], lo$psis_object, type = "mean", log_ratios = -ll[, obs])$value
  q  <- E_loo(pp[, obs], lo$psis_object, type = "quantile", probs = c(.05, .95), log_ratios = -ll[, obs])$value
  list(ep = ep, pp = pp, ll = ll, psis_obj = lo$psis_object, mu_loo = mu, q_loo = q, elpd = lo$estimates["elpd_loo", ], k_bad = sum(lo$diagnostics$pareto_k > 0.7),
       metrics = c(R2 = 1 - sum((y - mu)^2) / sum((y - mean(y))^2), RMSE = sqrt(mean((y - mu)^2)), RMSE_null = sd(y),
                   cover90 = mean(y >= q[1, ] & y <= q[2, ])), ep_sub = ep[sample(nd, 500), ], pp_sub = pp[sample(nd, 500), ])
}
# downstream_DOC_me: its own outcome, and dDOC derived
obsO <- which(!is.na(m6$DOC_out)); both <- which(!is.na(m6$DOC_out) & !is.na(m6$DOC_input))
O <- loo_block(fo, "DOCout", m6$DOC_out[obsO], obsO)
d_obs <- m6$DOC_out[both] - m6$DOC_input[both]
ib <- match(both, obsO); d_loo_O <- O$mu_loo[ib] - m6$DOC_input[both]
O$ddoc_metrics <- c(R2 = 1 - sum((d_obs - d_loo_O)^2) / sum((d_obs - mean(d_obs))^2), RMSE = sqrt(mean((d_obs - d_loo_O)^2)), RMSE_null = sd(d_obs))
ddO <- sweep(O$ep[, both], 2, m6$DOC_input[both])                                  # derived dDOC, draws x samplings
O$dd_med <- apply(ddO, 2, median); O$dd_lo <- apply(ddO, 2, quantile, .025); O$dd_hi <- apply(ddO, 2, quantile, .975)
O$avg_season <- sapply(c("summer", "winter"), function(s) rowMeans(ddO[, m6$season[both] == s]))
# LOO-adjusted 90% predictive interval for the derived dDOC (predicted DOC_out - observed
# upstream DOC), reweighted with the same PSIS weights as DOC_out's own leave-one-out.
dd_draws_obsO <- sweep(O$pp[, obsO], 2, ifelse(is.na(m6$DOC_input[obsO]), 0, m6$DOC_input[obsO]))  # placeholder 0 where DOC_input NA; those columns are dropped by [, ib] below
O$dd_q_loo <- E_loo(dd_draws_obsO, O$psis_obj, type = "quantile", probs = c(.05, .95), log_ratios = -O$ll[, obsO])$value[, ib]
O$ep <- NULL; O$pp <- NULL; O$ll <- NULL; O$psis_obj <- NULL
# delta_SCM_me: its own outcome is dDOC
obsD <- which(!is.na(m6$delta_DOC)); stopifnot(all(obsD == both))
D <- loo_block(fd, "deltaDOC", m6$delta_DOC[obsD], obsD)
D$dd_med <- apply(D$ep[, both], 2, median); D$avg_season <- sapply(c("summer", "winter"), function(s) rowMeans(D$ep[, both][, m6$season[both] == s]))
D$ep <- NULL
saveRDS(list(cmpD = cmpD, cmpO = cmpO, s_up = s_up, mu_up = mu_up, O = O, D = D, both = both, obsO = obsO, d_obs = d_obs,
             y_out = m6$DOC_out[obsO], up_both = m6$DOC_input[both], season_both = m6$season[both],
             rhat = c(delta = max(rhat(fd), na.rm = TRUE), down = max(rhat(fo), na.rm = TRUE))), here("results", "me_predictions.rds"))
print(round(O$metrics, 3)); print(round(O$ddoc_metrics, 3)); print(round(D$metrics, 3)); cat("done\n")
