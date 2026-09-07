# Posterior predictions, derived dDOC and PSIS-LOO for the downstream-DOC SCM
# (b1_downstream). Sourced by scm_outcome_choice.qmd and analysis_plan.qmd;
# results are cached in results/scm_downstream_pred.rds.
suppressMessages({library(brms); library(dplyr); library(loo); library(here)})
fit <- readRDS(here("results", "models_fit", "b1_downstream.rds"))
dat <- fit$data
n_chain <- 4; n_draw <- ndraws(fit)
ep <- posterior_epred(fit, resp = "DOCout")      # expected downstream DOC, draws x 360
pp <- posterior_predict(fit, resp = "DOCout")    # predicted new observation
ll <- log_lik(fit, resp = "DOCout")              # NA for the 22 rows without DOC_out
obs <- which(!is.na(dat$DOC_out))
both <- which(!is.na(dat$DOC_out) & !is.na(dat$DOC_input))
chain_id <- rep(seq_len(n_chain), each = n_draw / n_chain)
r_eff <- relative_eff(exp(ll[, obs]), chain_id = chain_id)
lo <- loo(ll[, obs], r_eff = r_eff, save_psis = TRUE)
psis <- lo$psis_object
# LOO-predictive mean and 90% interval per observed row
mu_loo <- E_loo(ep[, obs], psis, type = "mean", log_ratios = -ll[, obs])$value
q_loo  <- E_loo(pp[, obs], psis, type = "quantile", probs = c(.05, .95), log_ratios = -ll[, obs])$value
y <- dat$DOC_out[obs]
loo_metrics <- c(
  R2_loo   = 1 - sum((y - mu_loo)^2) / sum((y - mean(y))^2),
  RMSE_loo = sqrt(mean((y - mu_loo)^2)),
  cover90  = mean(y >= q_loo[1, ] & y <= q_loo[2, ]))
# derived dDOC: LOO-predictive downstream minus observed upstream
ib <- match(both, obs)
d_obs  <- dat$DOC_out[both] - dat$DOC_input[both]
d_loo  <- mu_loo[ib] - dat$DOC_input[both]
ddoc_metrics <- c(
  R2_loo   = 1 - sum((d_obs - d_loo)^2) / sum((d_obs - mean(d_obs))^2),
  RMSE_loo = sqrt(mean((d_obs - d_loo)^2)),
  RMSE_null = sqrt(mean((d_obs - mean(d_obs))^2)))
# in-sample Bayes R2 and posterior average dDOC by season (derived quantity)
# in-sample Bayes R2 computed by hand on the observed rows (bayes_R2() returns NA with mi() responses)
res <- sweep(-ep[, obs], 2, -dat$DOC_out[obs])            # y - mu per draw
r2_draw <- apply(ep[, obs], 1, var) / (apply(ep[, obs], 1, var) + apply(res, 1, var))
r2_in <- matrix(c(mean(r2_draw), sd(r2_draw), quantile(r2_draw, .025), quantile(r2_draw, .975)), nrow = 1,
                dimnames = list("R2DOCout", c("Estimate", "Est.Error", "Q2.5", "Q97.5")))
dd <- sweep(ep[, both], 2, dat$DOC_input[both])   # dDOC per draw and row
avg_season <- sapply(c("summer", "winter"), function(s) rowMeans(dd[, dat$season[both] == s]))
saveRDS(list(dat = dat, obs = obs, both = both,
             ep_sub = ep[sample(n_draw, 1000), ], pp_sub = pp[sample(n_draw, 1000), ],
             mu_loo = mu_loo, q_loo = q_loo, y = y,
             loo = lo, loo_metrics = loo_metrics, ddoc_metrics = ddoc_metrics,
             d_obs = d_obs, d_loo = d_loo, dd_sub = dd[sample(n_draw, 1000), ],
             avg_season = avg_season, r2_in = r2_in),
        here("results", "scm_downstream_pred.rds"))
print(lo); print(round(loo_metrics, 3)); print(round(ddoc_metrics, 3)); print(round(r2_in, 3))
