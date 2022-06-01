set.seed(42)

n_samples <- 100000

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

# Estimate impact of falls on health utility from Latimer et al (2013) Age and Ageing
# https://academic.oup.com/ageing/article/42/5/641/18607?login=true
noinj <- rep(0.02, 40)
minorinj <- rep(0.04, 31)
moderateinj <- rep(0.06, 18)
majorinj <- rep(0.11, 9)
fallers <- c(noinj, minorinj, moderateinj, majorinj)
falldec <- fitdist(fallers, "beta")

falls_params <- c(
  falls_params,
  list(fall_decrement=list(
    shape1=falldec$estimate['shape1'],
    shape2=falldec$estimate['shape2']
  ))
)

# ICU
icu_params <- list()

icu_params <- c(
  icu_params,
  list(opp_cost=436*1.03^(2022-2017))
)

# Estimate additional LOS of ICU readmission from Chen et al (1998) Crit Care Med
# https://pubmed.ncbi.nlm.nih.gov/9824076/
mean <- 7.8
sd <- 13.4
shape <- (mean^2)/(sd^2)
scale <- 1/(mean/(sd^2))
x <- rgamma(n_samples, shape = shape, scale = scale)
fit <- fitdist(x, "gamma", "mme")

icu_params <- c(
  icu_params,
  list(icu_readmit_los=list(
    shape=fit$estimate['shape'],
    rate=fit$estimate['rate'],
    multiplier=1.03^(2022-2021)
  ))
)

# Estimate impact of a day on the ward vs ICU from de Vos et al (2022), Value in Health
# https://www.sciencedirect.com/science/article/abs/pii/S1098301521017423

x <- rnorm(n_samples, 0.42, 0.083)
fit <- fitdist(x, "beta")

icu_params <- c(
  icu_params,
  list(ward_eff=list(
    shape1=fit$estimate['shape1'],
    shape2=fit$estimate['shape2'],
    multiplier=1/365
  ))
)

# Estimate ICU day cost from Hicks et al (2021), MJA
# https://www.mja.com.au/journal/2019/211/7/financial-cost-intensive-care-australia-multicentre-registry-study

x <- rnorm(n_samples, 4375, 1157)
fit <- fitdist(x[x>0], "gamma")

icu_params <- c(
  icu_params,
  list(icu_cost=list(
    shape=fit$estimate['shape'],
    rate=fit$estimate['rate'],
    multiplier=1.03^(2022-2021)
  ))
)

params <- list(global=global_params, falls=falls_params, icu=icu_params)

rm(list=setdiff(ls(), "params"))
