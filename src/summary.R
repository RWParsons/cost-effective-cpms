# summary functions

get_medians <- function(x, nboot=500){
  m <- median(x)

  boots <- map_dbl(
    1:nboot,
    function(v){
      s <- sample(x, size=length(v), replace=T)
      median(s)
    }
  )
  list(median=m, median_se=sd(boots))
}

summarize_sims <- function(x, prob=0.9, use_hdi=TRUE) {
  r <- function(x) round(x, digits=2)

  if(use_hdi){
    hdi <- hdi(x, ci=prob)
    res <- glue::glue(
      "{r(median(x))} [{scales::percent(prob)} HDI:{r(hdi$CI_low)}, {r(hdi$CI_high)}]"
    )
  } else {
    iqr <- quantile(x, probs=c(0.25,0.75))
    res <- glue::glue(
      "{r(median(x))} [IQR:{r(iqr[[1]])}, {r(iqr[[2]])}]"
    )
  }
  list(summary=res)
}

get_summary <- function(data, hdi_prob, use_hdi,
                        sample_size, n_sims, n_valid, sim_auc, event_rate) {
  df_summary <- data %>%
    pivot_longer(!n_sim, names_to="method", values_to="nmb")

  df_summary <- as.data.table(df_summary)[
    ,
    summarize_sims(x=nmb, prob=hdi_prob, use_hdi=use_hdi),
    by = list(method)
  ]

  df_summary$sample_size <- sample_size
  df_summary$n_sims <- n_sims
  df_summary$n_valid <- n_valid
  df_summary$sim_auc <- sim_auc
  df_summary$event_rate <- event_rate
  df_summary
}


# plotting functions

get_boxplot <- function(data, ordered_methods, subtitle) {
  data %>%
    pivot_longer(!n_sim, names_to="method", values_to="nmb") %>%
    mutate(method=factor(method, levels=ordered_methods)) %>%
    ggplot(aes(method, nmb)) +
    geom_boxplot() +
    geom_jitter(alpha=0.1, width=0.1) +
    theme_bw() +
    labs(
      y="", x="",
      subtitle=subtitle
    ) +
    scale_y_continuous(labels=scales::dollar_format()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

get_errorbar_plot <- function(data, ordered_methods, subtitle) {
  data <-
    data %>%
    pivot_longer(!n_sim, values_to="nmb", names_to="method")

  data <- as.data.table(data)[
    ,
    get_medians(x=nmb),
    by = list(method)
  ]

  data %>%
    mutate(method=factor(method, levels=ordered_methods)) %>%
    ggplot(aes(method, median)) +
    geom_point() +
    geom_errorbar(aes(ymin=median-median_se, ymax=median+median_se, width=0)) +
    theme_bw() +
    labs(
      y="", x="",
      subtitle=subtitle
    ) +
    scale_y_continuous(labels=scales::dollar_format()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}


plot_density_ridge = function(data, FUN=c("eti", "hdi"), ci=0.9, subtitle="", factor_levels=NULL, limit_y=FALSE, ...) {
  # taken from:
  # https://stackoverflow.com/questions/65269825/is-it-possible-to-recreate-the-functionality-of-bayesplots-mcmc-areas-plot-in

  # Determine whether to use eti or hdi function
  FUN = match.arg(FUN)
  FUN = match.fun(FUN)

  # Get kernel density estimate as a data frame
  dens = map_df(data, ~ {
    d = density(.x, na.rm=TRUE)
    tibble(x=d$x, y=d$y)
  }, .id="name")

  # Set relative width of median line
  e = diff(range(dens$x)) * 0.004

  # Get credible interval width and median
  cred.int = data %>%
    pivot_longer(cols=everything()) %>%
    group_by(name) %>%
    summarise(CI=list(FUN(value, ci=ci)),
              m=median(value, na.rm=TRUE)) %>%
    unnest_wider(CI)

  p_data <- left_join(dens, cred.int)

  if(!is.null(factor_levels)) {
    p_data$name <- factor(p_data$name, levels=factor_levels)
  }

  p <-
    p_data %>%
    ggplot(aes(y=name, x=x, height=y)) +
    geom_ridgeline(data= . %>% group_by(name) %>%
                     filter(between(x, CI_low, CI_high)),
                   fill=hcl(230,25,85), ...) +
    geom_ridgeline(data=. %>% group_by(name) %>%
                     filter(between(x, m - e, m + e)),
                   fill=hcl(240,30,60), ...) +
    geom_ridgeline(fill=NA, ...) +
    geom_ridgeline(fill=NA, aes(height=0), ...) +
    labs(y=NULL, x=NULL) +
    theme_bw() +
    coord_flip()

  if(limit_y){
    p <- p +
      scale_x_continuous(limits=c(0,1))
  }
  p +
    theme_bw() +
    labs(
      y="", x="",
      subtitle=subtitle
    ) +
    scale_x_continuous(labels=scales::dollar_format()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

# extract content from lists (made from parallel processing of simulations)
extract_plots <- function(l) {
  res <- list()
  for(i in 1:length(l)) {
    res <- c(res, list(l[[i]]$plot))
  }
  res
}

extract_summaries <- function(l) {
  res <- list()
  for(i in 1:length(l)) {
    res <- c(res, list(l[[i]]$summary))
  }
  res
}
