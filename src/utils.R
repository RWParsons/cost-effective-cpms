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

get_thresholds <- function(predicted, actual, pt_seq=seq(0.01, 0.99,0.01), costs){
  rocobj <- pROC::roc(as.factor(actual), predicted, direction="<", quiet=TRUE)
  auc <- pROC::auc(rocobj)
  pt_er <- pROC::coords(rocobj, "best", best.method="closest.topleft")$threshold
  pt_youden <- pROC::coords(rocobj, "best", best.method="youden")$threshold

  f <- function(pt){
    # new threshold selection methods are from here: https://www.hindawi.com/journals/cmmm/2017/3762651/
    cm <- get_confusion(d=data.frame(predicted=predicted, actual=actual), pt=pt)
    data.frame(
      pt=pt,
      cost_effective=cm$TN*costs["TN"] + cm$TP*costs["TP"] + cm$FN*costs["FN"] + cm$FP*costs["FP"],
      cz=cm$Se*cm$Sp,
      iu=abs(cm$Se - auc) + abs(cm$Sp - auc)
    )
  }

  screen_df <- map_dfr(pt_seq, f)

  pt_cost_effective <- mean(screen_df$pt[screen_df$cost_effective==min(screen_df$cost_effective)])
  pt_cz <- mean(screen_df$pt[screen_df$cz==max(screen_df$cz)])
  pt_iu <- mean(screen_df$pt[screen_df$iu==min(screen_df$iu)])

  list(pt_er=pt_er, pt_youden=pt_youden, pt_cost_effective=pt_cost_effective, pt_cz=pt_cz, pt_iu=pt_iu)
}

get_confusion <- function(d, pt=0.2){
  TN <- sum(d$predicted < pt & d$actual==0)
  FN <- sum(d$predicted < pt & d$actual==1)
  TP <- sum(d$predicted > pt & d$actual==1)
  FP <- sum(d$predicted > pt & d$actual==0)

  Se <- TP/(TP+FN)
  Sp <- TN/(FP+TN)

  list(TN=TN, FN=FN, TP=TP, FP=FP, Se=Se, Sp=Sp)
}
