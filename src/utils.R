get_sample <- function(auc, n_samples, prevalence, scale_to_d=F){
  # https://stats.stackexchange.com/questions/422926/generate-synthetic-data-given-auc
  # http://dx.doi.org/10.5093/ejpalc2018a5
  t <- sqrt(log(1/(1-auc)**2))
  z <- t-((2.515517 + 0.802853*t + 0.0103328*t**2) /
            (1 + 1.432788*t + 0.189269*t**2 + 0.001308*t**3))
  d <- z*sqrt(2)

  n_pos <- sum(sample(c(0,1), n_samples, replace=TRUE, prob=c(1-prevalence, prevalence)))
  n_neg <- n_samples - n_pos

  x <- c(rnorm(n_neg, mean=0), rnorm(n_pos, mean=d))
  y <- c(rep(0, n_neg), rep(1, n_pos))

  if(scale_to_d){
    x <- x/d
  }
  return(data.frame(predicted=x, actual=y))
}

# Beta distribution-related functions for simulating clinical prediction models

# get_alpha() and get_beta() calculate one given the other and the the prevalence (p)
# this exploits the formula for the expectation of the beta distribution:
# E[X] = a/(a+b)
# this assumes that the expectation of predicted probabilities is equal to the prevalence of the event
get_alpha <- function(beta, p){
  -(beta*p)/(p-1)
}

get_beta <- function(alpha, p){
  (alpha-alpha*p)/p
}

get_auc <- function(predicted, actual){
  AUC::auc(AUC::roc(predicted, as.factor(actual)))
}

get_beta_preds <- function(alpha=NULL, beta=NULL, p=NULL, n, get_what=c("preds", "auc", "params")){
  if(is.null(alpha)){
    alpha <- get_alpha(beta=beta, p=p)
  }
  if(is.null(beta)){
    beta <- get_beta(alpha=alpha, p=p)
  }

  predicted_probs <- rbeta(n=n, shape1=alpha, shape2=beta)
  f <- function(x) sample(c(0, 1), 1, prob=c(1-x, x))
  predicted_classes <- purrr::map_dbl(predicted_probs, f)

  res <- list()
  if("params" %in% get_what){
    res <- list(alpha=alpha, beta=beta)
  }
  if("preds" %in% get_what){
    if(length(res)==0){
      res <- list(preds=data.frame(predicted=predicted_probs, actual=predicted_classes))
    }else{
      res <- append(res, list(preds=data.frame(predicted=predicted_probs, actual=predicted_classes)))
    }
  }
  if("auc" %in% get_what){
    if(length(res)==0){
      res <- list(auc=get_auc(predicted=predicted_probs, actual=predicted_classes))
    }else{
      res <- append(res, list(auc=get_auc(predicted=predicted_probs, actual=predicted_classes)))
    }
  }
  res
}

# given some predictions corresponding labels, a probability threshold and a vector containing costs, calculate the total cost
classify_samples <- function(predicted, actual, pt, costs){
  # predicted: vector of predicted probabilities
  # actual: binary label to the event
  # pt: probability threshold used to classify predicted probabilities into classes
  # costs:  named vector containing costs for each possible correct or incorrect classification (2x2)
  #         for example: costs <- c("TN"=0, "FN"=100, "TP"=80, "FP"=5)

  d <- cbind(predicted, actual, NA)
  colnames(d) <- c("predicted", "actual", "nmb")

  d[d[,"predicted"] < pt & d[,"actual"]==0, "nmb"] <- costs["TN"]

  d[d[,"predicted"] < pt & d[,"actual"]==1, "nmb"] <- costs["FN"]
  d[d[,"predicted"] > pt & d[,"actual"]==1, "nmb"] <- costs["TP"]
  d[d[,"predicted"] > pt & d[,"actual"]==0, "nmb"] <- costs["FP"]

  mean(d[,"nmb"])
}

get_smooth_max <- function(x, y){
  smooth <- do.call(supsmu, list(x, y))
  max.idx <- which.max(smooth$y)
  smooth$x[max.idx]
}


get_thresholds <- function(predicted, actual, costs){
  pt_er <- cutpointr(
    x=predicted, class=actual, method=minimize_metric, metric=roc01,
    silent=TRUE
  )[['optimal_cutpoint']]
  if(length(pt_er) > 1) {
    pt_er <- median(pt_er)
  }

  pt_youden <- cutpointr(
    x=predicted, class=actual, method=maximize_metric, metric=youden,
    silent=TRUE
  )[['optimal_cutpoint']]
  if (pt_youden > 1) {
    pt_youden <- 1
  }

  pt_cost_effective <- cutpointr(
    x=predicted, class=actual, method=maximize_metric, metric=fx_total_nmb,
    utility_tp=costs["TP"], utility_tn=costs["TN"],
    cost_fp=costs["FP"], cost_fn=costs["FN"],
    silent=TRUE
  )[['optimal_cutpoint']]
  if (pt_cost_effective > 1) {
    pt_cost_effective <- 1
  }

  pt_cz <- cutpointr(
    x=predicted, class=actual, method=maximize_metric, metric=fx_prod_sn_sp,
    silent=TRUE
  )[['optimal_cutpoint']]
  if (pt_cz > 1) {
    pt_cz <- 1
  }

  pt_iu <- cutpointr(
    x=predicted, class=actual, method=minimize_metric, metric=roc_iu,
    silent=TRUE
  )[['optimal_cutpoint']]
  if (pt_iu > 1) {
    pt_iu <- 1
  }

  costs_pos <- -costs

  pt_cost_minimising <-
    (costs_pos["FP"] - costs_pos["TN"]) /
    (costs_pos["FP"] + costs_pos["FN"] - costs_pos["TP"] - costs_pos["TN"])
  names(pt_cost_minimising) <- NULL

  list(
    treat_all=0,
    treat_none=1,
    cost_effective=pt_cost_effective,
    er=pt_er,
    youden=pt_youden,
    cz=pt_cz,
    iu=pt_iu,
    cost_minimising=pt_cost_minimising
  )
}

get_confusion <- function(predicted, actual, pt){
  TN <- sum(predicted < pt & actual==0)
  FN <- sum(predicted < pt & actual==1)
  TP <- sum(predicted > pt & actual==1)
  FP <- sum(predicted > pt & actual==0)

  Se <- TP/(TP+FN)
  Sp <- TN/(FP+TN)

  list(TN=TN, FN=FN, TP=TP, FP=FP, Se=Se, Sp=Sp)
}
