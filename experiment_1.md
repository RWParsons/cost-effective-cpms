Experiment 1
================
15 February, 2022

Question: What are the differences in NMB between models where the
Probability threshold was based on the currently available methods
versus costs-based selection. (Hospital falls as a use case.)

1.  Define costs of a TP, TN, FP, FN of falls classification (option to
    move this into the loop where costs are sampled from a distributions
    to account for uncertainty in their estimates in the literature)
    -   FP have cost of applying intervention
    -   FN have cost of patient fall
    -   TP have cost of intervention + cost of fall\*(1-effectiveness of
        intervention on rate of falls)
    -   TN are cost $0
2.  Select appropriate AUC (\~0.75?) and prevalence (\~3% ) for
    comparable clinical prediction model for falls.
3.  For sample sizes (N) in \[100, 500, 1000\]: (repeat 500 times at
    each sample size)
    -   Get training data by sampling observed predictor values and
        outcome by transforming AUC into Cohens’ D and sampling from two
        normal distributions, the first (negative events) with mean=0
        and the second (positive events) with mean=Cohens’D. (Both with
        sd=1.)
    -   Fit a logistic regression model using this sampled data.
    -   Fit predicted probabilities to the training data and use these
        to obtain probability thresholds using each method.
    -   Get validation data using the same approach but with n=1000.
    -   Use the previously fit model to estimate probabilities for
        validation data.
    -   Evaluate the thresholds selected using the training data on the
        validation data, in terms of mean cost per patient.
4.  Measure differences in NMB on validation sample dependent on use of
    currently available methods and cost-based approach to determine
    threshold.
5.  Observe whether this relationship is dependent on the sample size
    taken

### Define costs

``` r
treatment_cost <-250
treatment_effect <- 0.1
wtp <- 28000
QALYs_event <- 0.1

cost_vector <- c(
  "TN"=0, 
  "FN"=QALYs_event*wtp, 
  "TP"=QALYs_event*(1-treatment_effect)*wtp, 
  "FP"=treatment_cost
)

f_plot <- function(beta, rate){
  alpha <- get_alpha(beta, rate)
  print(alpha)
  data.frame(sample=rbeta(n=10000, alpha, beta)) %>%
    ggplot(aes(sample)) + geom_density()
}
# f_plot(1000, 0.1)

get_costs <- function(){
  treatment_effect <- rbeta(1, 111.111, 1000)
  QALYs_event <- rbeta(1, 10/3, 30)
  treatment_cost <- rgamma(1, 250)
  wtp <- 28000
  
  c(
    "TN"=0, 
    "FN"=QALYs_event*wtp, 
    "TP"=QALYs_event*(1-treatment_effect)*wtp, 
    "FP"=treatment_cost
  )
}
get_costs()
```

    ##        TN        FN        TP        FP 
    ##    0.0000 3504.4305 3173.8869  259.5966

### Add helper functions for determining classifying predictions and estimating probability thresholds.

``` r
get_confusion <- function(d, pt=0.2){
  TN <- sum(d$predicted < pt & d$actual==0)
  FN <- sum(d$predicted < pt & d$actual==1)
  TP <- sum(d$predicted > pt & d$actual==1)
  FP <- sum(d$predicted > pt & d$actual==0)
  
  Se <- TP/(TP+FN)
  Sp <- TN/(FP+TN)
  
  list(TN=TN, FN=FN, TP=TP, FP=FP, Se=Se, Sp=Sp)
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
```

### Run simulation

``` r
# source("src/utils.R")
sample_sizes <- c(100, 500, 1000)
n_sims <- 100
n_valid <- 1000
sim_auc <- 0.65
event_rate <- 0.02

df_result <- data.frame(
  train_sample_size=c(), n_sim=c(), youden_cost=c(), cost_effective_cost=c(),
  cz_cost=c(), iu_cost=c(), er_cost=c()
)


for(n in sample_sizes){
  i <- 0
  while(i < n_sims){
    train_sample <- get_sample(auc=sim_auc, n_samples=n, prevalence=event_rate)
    valid_sample <- get_sample(auc=sim_auc, n_samples=n_valid, prevalence=event_rate)
    if(length(unique(train_sample$actual))!=2 | length(unique(valid_sample$actual))!=2){
      next
    }
    i <- i + 1
    model <- glm(actual~predicted, data=train_sample, family=binomial())
    train_sample$predicted <- predict(model, type="response")

    
    valid_sample$predicted <- predict(model, type="response", newdata=valid_sample)
    
    cost_vector <- get_costs()
    
    thresholds <- get_thresholds(
      predicted=train_sample$predicted, 
      actual=train_sample$actual,
      costs=cost_vector
    )
    
    youden_mean_cost <- classify_samples(
      predicted=valid_sample$predicted,
      actual=valid_sample$actual,
      pt=thresholds$pt_youden,
      costs=cost_vector
    )
    
    cost_vector <- get_costs()
    
    cost_threshold <- function(pt){
      classify_samples(
        predicted=valid_sample$predicted,
        actual=valid_sample$actual,
        pt=pt,
        costs=cost_vector
      )
    }
    
    cost_effective_mean_cost <- classify_samples(
      predicted=valid_sample$predicted,
      actual=valid_sample$actual,
      pt=thresholds$pt_cost_effective,
      costs=cost_vector
    )
    
    cost_effective_mean_cost <- classify_samples(
      predicted=valid_sample$predicted,
      actual=valid_sample$actual,
      pt=thresholds$pt_cost_effective,
      costs=cost_vector
    )
    
    df_result <- rbind(
      df_result, 
      data.frame(
        train_sample_size=n, 
        n_sim=i, 
        youden_cost=cost_threshold(thresholds$pt_youden), 
        cost_effective_cost=cost_threshold(thresholds$pt_cost_effective),
        cz_cost=cost_threshold(thresholds$pt_cz), 
        er_cost=cost_threshold(thresholds$pt_er), 
        iu_cost=cost_threshold(thresholds$pt_iu)
      )
    )
  }
}

df_result %>%
  select(-n_sim) %>%
  pivot_longer(!train_sample_size, names_to="method", values_to="cost") %>%
  mutate(method=str_remove(method, "_cost")) %>%
  ggplot(aes(method, cost)) +
  geom_boxplot() +
  geom_jitter(alpha=0.1, width=0.1) +
  facet_wrap(~train_sample_size) +
  theme_bw() +
  labs(y="Mean cost per patient") +
  scale_y_continuous(labels=scales::dollar_format()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](experiment_1_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
