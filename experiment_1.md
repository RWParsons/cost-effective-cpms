Experiment 1
================
18 February, 2022

Question: What are the differences in NMB between models where the
Probability threshold was based on the currently available methods
versus costs-based selection. (Hospital falls as a use case.)

1.  Define costs of a TP, TN, FP, FN of falls classification (option to
    move this into the loop where costs are sampled from a distributions
    to account for uncertainty in their estimates in the literature)
      - FP have cost of applying intervention
      - FN have cost of patient fall
      - TP have cost of intervention + cost of fall\*(1-effectiveness of
        intervention on rate of falls)
      - TN are cost $0
2.  Select appropriate AUC (\~0.75?) and prevalence (\~3% ) for
    comparable clinical prediction model for falls.
3.  For sample sizes (N) in \[100, 500, 1000\]: (repeat 500 times at
    each sample size)
      - Get training data by sampling observed predictor values and
        outcome by transforming AUC into Cohens’ D and sampling from two
        normal distributions, the first (negative events) with mean=0
        and the second (positive events) with mean=Cohens’D. (Both with
        sd=1.)
      - Fit a logistic regression model using this sampled data.
      - Fit predicted probabilities to the training data and use these
        to obtain probability thresholds using each method.
      - Get validation data using the same approach but with n=1000.
      - Use the previously fit model to estimate probabilities for
        validation data.
      - Evaluate the thresholds selected using the training data on the
        validation data, in terms of mean cost per patient.
4.  Measure differences in NMB on validation sample dependent on use of
    currently available methods and cost-based approach to determine
    threshold.
5.  Observe whether this relationship is dependent on the sample size
    taken
6.  ???

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
    "TP"=QALYs_event*(1-treatment_effect)*wtp + treatment_cost, 
    "FP"=treatment_cost
  )
}
get_costs()
```

    ##        TN        FN        TP        FP 
    ##    0.0000 3504.4305 3433.4835  259.5966

``` r
get_nmb <- function(){
  treatment_effect <- rbeta(1, 111.111, 1000) #QALY increment from treatment
  QALYs_event <- rbeta(1, 10/3, 30)*-1 #QALY decrement from event
  treatment_cost <- rgamma(1, 250)
  wtp <- 28000
  
  c(
    "TN"=0,
    "FN"=QALYs_event*wtp,
    "TP"=QALYs_event*(1-treatment_effect)*wtp - treatment_cost,
    "FP"=treatment_cost*-1
  )
}
```

### Run simulation

``` r
do_simulation <- function(sample_size, n_sims, n_valid, sim_auc, event_rate, fx_costs, get_what=c("data", "plot")){
  df_result <- data.frame()
  df_thresholds <- data.frame()
  
  i <- 0
  while(i < n_sims){
    train_sample <- get_sample(auc=sim_auc, n_samples=sample_size, prevalence=event_rate)
    valid_sample <- get_sample(auc=sim_auc, n_samples=n_valid, prevalence=event_rate)
    if(length(unique(train_sample$actual))!=2 | length(unique(valid_sample$actual))!=2){
      next
    }
    i <- i + 1
    model <- glm(actual~predicted, data=train_sample, family=binomial())
    train_sample$predicted <- predict(model, type="response")

    
    valid_sample$predicted <- predict(model, type="response", newdata=valid_sample)
    
    cost_vector <- fx_costs()
    
    thresholds <- get_thresholds(
      predicted=train_sample$predicted, 
      actual=train_sample$actual,
      costs=cost_vector
    )
    
    cost_vector <- fx_costs()
    
    cost_threshold <- function(pt){
      classify_samples(
        predicted=valid_sample$predicted,
        actual=valid_sample$actual,
        pt=pt,
        costs=cost_vector
      )
    }
    
    df_result <- rbind(
      df_result, 
      data.frame(
        n_sim=i, 
        treat_all=cost_threshold(0),
        treat_none=cost_threshold(1),
        cost_effective_cost=cost_threshold(thresholds$pt_cost_effective),
        youden_cost=cost_threshold(thresholds$pt_youden), 
        cz_cost=cost_threshold(thresholds$pt_cz), 
        iu_cost=cost_threshold(thresholds$pt_iu),
        er_cost=cost_threshold(thresholds$pt_er)
        
      )
    )
    
    df_thresholds <- rbind(
      df_thresholds, 
      data.frame(
        n_sim=i, 
        treat_all=0,
        treat_none=1,
        cost_effective_pt=thresholds$pt_cost_effective,
        youden_pt=thresholds$pt_youden, 
        cz_pt=thresholds$pt_cz,
        iu_pt=thresholds$pt_iu,
        er_pt=thresholds$pt_er
      )
    )
  }
  
  pts <- round(colMeans(df_thresholds)[-1], 3)
  df_plot <- df_result
  names(df_plot)[-1] <- method_levels <- paste0(names(df_plot)[-1], "(", pts, ")")
  
  p <- 
    df_plot %>%
    pivot_longer(!n_sim, names_to="method", values_to="cost") %>%
    mutate(method=str_remove(method, "_cost"),
           method=factor(method, levels=str_remove(method_levels, "_cost"))) %>%
    
    ggplot(aes(method, cost)) +
    geom_boxplot() +
    geom_jitter(alpha=0.1, width=0.1) +
    theme_bw() +
    labs(
      # y="Mean cost per patient",
      # x="Probability threshold selection method",
      y="", x="",
      subtitle=glue::glue("auc: {sim_auc}; event_rate: {event_rate}")
    ) +
    scale_y_continuous(labels=scales::dollar_format()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  if(length(get_what)==1){
    if(get_what=="data"){return(df_result)}
    if(get_what=="plot"){return(p)}
  }else{
    return(list(data=df_result, plot=p))
  }
}

# do_simulation(sample_size=500, n_sims=20, n_valid=1000, sim_auc=0.65, event_rate=0.03, fx_costs=get_costs, get_what=c("data", "plot"))

# p_list <- lapply(c(0.03, 0.1, 0.5, 0.9, 0.97), function(x)do_simulation(sample_size=500, n_sims=100, n_valid=1000, sim_auc=0.65, event_rate=x, fx_costs=get_costs, get_what="plot"))
# p_list
```

## search through a grid of combinations of AUC and event rates to see how this influences the differences between probability threshold methods. The same costs were used in all simulations (distributions at top of document) and are resampled separately for training and validation.

#### In the plot below, the columns, from left to right, have increasing AUC. The rows, from top to bottom, have increasing event rates.

### My hot takes:

#### Probability threshold selection method becomes increasingly important as the AUC of the model and the event rate reduces.

#### I think that this is because, for models with very high discrimination, they’re able to correctly classify a larger proportion of samples, and there are fewer which are classified differently based on selection method. The difference is greater for smaller event rates because a false negative is the most costly classification, and only the cost-based method is “aware” of this. This is also why, when the event rate is very high, there is not much of a difference between methods (cost-based method is focused on correctly classifying the majority class but so are the other methods).

``` r
library(parallel)
n_cluster <- detectCores()
cl <- makeCluster(n_cluster)
cl <- parallelly::autoStopCluster(cl)

g <- expand.grid(
  sim_auc=c(0.65, 0.75, 0.85, 0.95),
  event_rate=c(0.01, 0.03, 0.1, 0.3, 0.7, 0.9)
)

clusterExport(cl, {
  c("do_simulation", "g", "get_costs")
})

invisible(clusterEvalQ(cl, {
  library(tidyverse)
  source("src/utils.R")
}))


ll <- parallel::parLapply(
  cl,
  1:nrow(g),
  function(i) do_simulation(
    sample_size=500, n_sims=100, n_valid=1000,
    sim_auc=g$sim_auc[i], event_rate=g$event_rate[i],
    fx_costs=get_costs, get_what="plot"
  )
)

cowplot::plot_grid(plotlist=ll, ncol=length(unique(g$sim_auc)))
```

![](experiment_1_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
