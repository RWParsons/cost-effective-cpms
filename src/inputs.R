set.seed(42)
library(tidyverse)
library(fitdistrplus)
n_samples <- 1000000

params <- list()

global_params <- list(
  WTP=28033
)

# Falls
falls_params <- list()

falls_params <- c(
  falls_params,
  list(treatment_cost=294*1.03^(2022-2013))
)

#### costs of falls from Morello et al. MJA (https://www.mja.com.au/journal/2015/203/9/extra-resource-burden-hospital-falls-cost-falls-study)
# the confidence interval from the study were symmetric around the estimate.
# The code below is used to estimate which gamma distribution parameters would represent costs with
# the same mean and standard error.
mean_cost <- 6669
ul <- 9450
desired_se <- ((ul-mean_cost)/1.96)

x <- rnorm(n=n_samples, mean_cost, desired_se)
x <- x[x>0]
fit <- fitdist(x, "gamma", method="mme")

falls_params <- c(
  falls_params,
  list(falls_cost=list(
    shape=fit$estimate['shape'],
    rate=fit$estimate['rate'],
    multiplier=(1.03^(2022-2015))
  ))
)

#### estimate effectiveness of intervention from Haines et al. (2010) Archives of Internal Medicine (https://sci-hub.hkvisa.net/10.1001/archinternmed.2010.444)
mean_eff <- 0.43
ci95 <- c(0.24, 0.78)

log_hr <- log(mean_eff)
log_hr_se <- (log(ci95[2]) - log(mean_eff))/1.96

falls_params <- c(
  falls_params,
  list(treatment_log_hazard=list(
    mean=log_hr,
    sd=log_hr_se
  ))
)


# ICU params
icu_params <- list()


params <- list(global=global_params, falls=falls_params, icu=icu_params)
