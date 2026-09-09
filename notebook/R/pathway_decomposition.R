# Pathway decomposition of the local modification predicted by the downstream-DOC
# SCM (b1_downstream): for every sampling, season and posterior draw, how much of
# (predicted downstream DOC - upstream DOC) comes from each modelled pathway.
# Sourced by analysis_plan.qmd; summaries cached in results/pathway_decomposition.rds.
suppressMessages({library(brms); library(dplyr); library(tidyr); library(here); library(posterior)})
fit <- readRDS(here("results", "models_fit", "b1_downstream.rds"))
dat <- fit$data
set.seed(7); nd <- 1000; keep <- sort(sample(ndraws(fit), nd))
dr <- as_draws_df(fit) |> subset_draws(draw = keep)
g <- function(pat) as.matrix(dr[, grep(pat, names(dr), value = TRUE), drop = FALSE])
# --- coefficients (draws x 1) ---
a   <- dr$b_DOCout_Intercept
b_in <- dr$bsp_DOCout_miDOC_input
b   <- list(phytoplankton = dr$bsp_DOCout_miplankton_abun_z, macrophytes = dr$bsp_DOCout_mimacrophy_abun_z,
            litter = dr$bsp_DOCout_micover_litter_z, hydrology = dr$bsp_DOCout_miwater_res_time_z)
# --- mediator values per row and draw: observed, or the draw's imputed value ---
fill <- function(var, resp) {
  X <- matrix(dat[[var]], nrow = nd, ncol = nrow(dat), byrow = TRUE)
  mis <- which(is.na(dat[[var]]))
  if (length(mis)) { Y <- g(paste0("^Ymi_", resp, "\\[")); stopifnot(ncol(Y) == length(mis)); X[, mis] <- Y }
  X }
M <- list(phytoplankton = fill("plankton_abun_z", "planktonabunz"), macrophytes = fill("macrophy_abun_z", "macrophyabunz"),
          litter = fill("cover_litter_z", "coverlitterz"), hydrology = fill("water_res_time_z", "waterrestimez"))
UP <- fill("DOC_input", "DOCinput")
# --- site effects per row ---
rs <- g("^r_site__DOCout\\["); sites <- sub("^r_site__DOCout\\[(.*),Intercept\\]$", "\\1", colnames(rs))
U <- rs[, match(as.character(dat$site), sites)]
# --- components of the local modification (draws x rows), mg/L ---
comp <- list(inheritance = sweep(UP, 1, b_in - 1, `*`))          # (beta_inh - 1) * DOC_up
for (k in names(b)) comp[[k]] <- sweep(M[[k]], 1, b[[k]], `*`)  # beta_k * M_k
comp$unresolved <- matrix(a, nd, nrow(dat)) + U                   # intercept + site effect
L <- Reduce(`+`, comp)                                            # predicted local modification = mu - DOC_up
# --- classes per sampling from the derived local modification ---
p_pos <- colMeans(L > 0)
cls <- ifelse(p_pos > 0.90, "source", ifelse(p_pos < 0.10, "sink", "no-evidence"))
measured <- !is.na(dat$plankton_abun_z)
rows <- tibble(row = seq_len(nrow(dat)), site = as.character(dat$site), season = as.character(dat$season),
               beaver_territory = dat$beaver_territory, measured = measured, class = cls,
               L_med = apply(L, 2, median), obs_ddoc = dat$DOC_out - dat$DOC_input)
# --- posterior of the mean contribution by class (all rows, and measured rows only) ---
by_class <- function(idx) bind_rows(lapply(names(comp), function(k) {
  X <- comp[[k]]
  bind_rows(lapply(c("source", "sink", "no-evidence"), function(cl) {
    j <- idx[cls[idx] == cl]; if (length(j) < 2) return(NULL)
    v <- rowMeans(X[, j, drop = FALSE]); tibble(component = k, class = cl, n = length(j), med = median(v), lo = quantile(v, .025), hi = quantile(v, .975)) })) }))
summ_all <- by_class(seq_len(nrow(dat))); summ_meas <- by_class(which(measured))
# source-minus-sink difference per component (posterior)
diff_all <- bind_rows(lapply(names(comp), function(k) { X <- comp[[k]]
  js <- which(cls == "source"); jk <- which(cls == "sink")
  v <- rowMeans(X[, js, drop = FALSE]) - rowMeans(X[, jk, drop = FALSE])
  tibble(component = k, med = median(v), lo = quantile(v, .025), hi = quantile(v, .975), p_pos = mean(v > 0)) }))
# --- seasonal change: summer minus winter per stream, per component, and B ---
pairs <- rows |> select(row, site, beaver_territory, season) |> pivot_wider(names_from = season, values_from = row) |>
  filter(!is.na(summer), !is.na(winter))
Bm <- L[, pairs$summer] - L[, pairs$winter]                       # B per draw and stream (from the SCM)
pB <- colMeans(Bm > 0); clsB <- ifelse(pB > 0.90, "B > 0 (summer > winter)", ifelse(pB < 0.10, "B < 0 (winter > summer)", "no-evidence"))
seas <- bind_rows(lapply(names(comp), function(k) { Dk <- comp[[k]][, pairs$summer] - comp[[k]][, pairs$winter]
  bind_rows(lapply(unique(clsB), function(cl) { j <- which(clsB == cl); if (length(j) < 2) return(NULL)
    v <- rowMeans(Dk[, j, drop = FALSE]); tibble(component = k, classB = cl, n = length(j), med = median(v), lo = quantile(v, .025), hi = quantile(v, .975)) })) }))
# share of B explained by each component: posterior of the slope of component change on B across streams
share <- bind_rows(lapply(names(comp), function(k) { Dk <- comp[[k]][, pairs$summer] - comp[[k]][, pairs$winter]
  v <- sapply(seq_len(nd), function(d) { x <- Bm[d, ]; y <- Dk[d, ]; sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2) })
  tibble(component = k, med = median(v), lo = quantile(v, .025), hi = quantile(v, .975)) }))
# variance share of the local modification across samplings, per draw
vshare <- bind_rows(lapply(names(comp), function(k) { v <- sapply(seq_len(nd), function(d) cov(comp[[k]][d, ], L[d, ]) / var(L[d, ]))
  tibble(component = k, med = median(v), lo = quantile(v, .025), hi = quantile(v, .975)) }))
saveRDS(list(rows = rows, summ_all = summ_all, summ_meas = summ_meas, diff_all = diff_all, seas = seas, share = share,
             vshare = vshare, n_class = table(cls), n_classB = table(clsB), n_pairs = nrow(pairs), nd = nd),
        here("results", "pathway_decomposition.rds"))
print(table(cls)); print(vshare); print(diff_all); cat("done\n")
