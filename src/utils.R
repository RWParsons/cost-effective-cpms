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

  d <- data.frame(
    predicted=predicted,
    actual=actual
  )
  d$cost <- NA
  d$cost[d$predicted < pt & d$actual==0] <- costs["TN"]
  d$cost[d$predicted < pt & d$actual==1] <- costs["FN"]
  d$cost[d$predicted > pt & d$actual==1] <- costs["TP"]
  d$cost[d$predicted > pt & d$actual==0] <- costs["FP"]

  mean(d$cost)
}


get_thresholds <- function(predicted, actual, costs, pt_seq=seq(0.01, 0.99,0.01)){
  df_pt_costs <- data.frame(pt=pt_seq)
  df_pt_costs$mean_cost <- map_dbl(
    df_pt_costs$pt,
    function(x)classify_samples(predicted=predicted, actual=actual, pt=x, costs=costs)
  )

  cost_effective_pt <- slice(arrange(df_pt_costs, mean_cost),1)$pt

  rocobj <- pROC::roc(as.factor(actual), predicted, direction="<", quiet=TRUE)
  youden_pt <- pROC::coords(rocobj, "best")$threshold
  res <- list(youden=youden_pt, cost_effective=cost_effective_pt)
  return(res)
}
