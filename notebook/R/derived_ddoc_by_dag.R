# Does the derived dDOC depend on the SCM specification? Recompute it from each
# of the four fitted DAGs (section 9 of analysis_plan.qmd) and compare.
suppressMessages({library(brms); library(dplyr); library(here)})
m6 <- read.csv(here("data/processed/m6.csv"))
out <- list()
for (nm in c("dag0","dag1","dag2","dag3")) {
  fit <- readRDS(here("results","models_fit", paste0("scm_", nm, ".rds")))
  dat <- fit$data; stopifnot(nrow(dat) == nrow(m6))
  both <- which(!is.na(dat$DOC_out) & !is.na(dat$DOC_input))
  set.seed(7); ep <- posterior_epred(fit, resp = "DOCout", ndraws = 1000)
  dd <- sweep(ep[, both], 2, dat$DOC_input[both])          # derived dDOC per draw and sampling
  out[[nm]] <- list(med = apply(dd, 2, median), lo = apply(dd, 2, quantile, .025), hi = apply(dd, 2, quantile, .975),
                    avg = sapply(c("summer","winter"), function(s) rowMeans(dd[, m6$season[both] == s])),
                    obs = dat$DOC_out[both] - dat$DOC_input[both], both = both)
  cat(nm, "done\n")
}
saveRDS(out, here("results", "derived_ddoc_by_dag.rds"))
