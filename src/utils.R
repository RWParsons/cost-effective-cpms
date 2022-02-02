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

get_beta_preds <- function(alpha=NULL, beta=NULL, p=NULL, n, return_preds=FALSE){
  if(is.null(alpha)){
    alpha <- get_alpha(beta=beta, p=p)
  }
  if(is.null(beta)){
    beta <- get_beta(alpha=alpha, p=p)
  }

  predicted_probs <- rbeta(n=n, shape1=alpha, shape2=beta)
  f <- function(x) sample(c(0, 1), 1, prob=c(1-x, x))
  predicted_classes <- purrr::map_dbl(predicted_probs, f)

  if(return_preds){
    return(data.frame(predicted=predicted_probs, actual=predicted_classes))
  }
  return(AUC::auc(AUC::roc(predicted_probs, as.factor(predicted_classes))))
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
