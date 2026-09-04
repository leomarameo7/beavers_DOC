# Latent (measurement-error corrected) DOC concentrations per reach, used for S, B, G and R_s.
# One hierarchical model per measured concentration, on the log scale (positive, right-skewed),
# with the 0.2 mg/L instrument error transferred by the delta method. Output: results/latent_R.rds
suppressMessages({library(brms);library(tidyverse);library(here)})
m <- read.csv(here("data","processed","m12.csv"))
w <- m|>filter(season=="winter")|>select(site,beaver_territory,up_w=DOC_input,dn_w=DOC_out,Afp=A_floodplain_km2,Abv=A_beaver_km2)
s <- m|>filter(season=="summer")|>select(site,beaver_territory,up_s=DOC_input,dn_s=DOC_out)
d <- inner_join(w,s,by=c("site","beaver_territory")) |> filter(complete.cases(up_w,dn_w,up_s,dn_s)) |> mutate(reach=factor(row_number()))
pri <- c(prior(normal(1,1.5), class="Intercept"), prior(student_t(3,0,1), class="sd"))
fitc <- function(v){ dd <- d |> mutate(y=log(pmax(.data[[v]],0.05)), se_log = 0.2/pmax(.data[[v]],0.25))
  brm(bf(y | se(se_log, sigma=FALSE) ~ 1 + (1|reach)), data=dd, prior=pri, chains=4, cores=4, iter=6000, warmup=2000,
      seed=7, backend="cmdstanr", control=list(adapt_delta=0.99, max_treedepth=12),
      file=here("results","models_fit",paste0("latent_",v)), refresh=0) }
lat <- function(f){x<-as_draws_df(f); exp(as.matrix(x[,grep("^r_reach\\[",names(x))]) + x$b_Intercept)}
L <- setNames(lapply(c("up_w","dn_w","up_s","dn_s"), function(v) lat(fitc(v))), c("up_w","dn_w","up_s","dn_s"))
S <- L$up_s - L$up_w; dw <- L$dn_w - L$up_w; ds <- L$dn_s - L$up_s
saveRDS(list(S=S,dw=dw,ds=ds,R_s=ds/abs(S),R_w=dw/abs(S),up_s=L$up_s,up_w=L$up_w,d=d), here("results","latent_R.rds"))
