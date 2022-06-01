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
      "{r(median(x))} [{r(hdi$CI_low)}, {r(hdi$CI_high)}]"
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


density_plot <- function(data, ci=0.95, limit_y=FALSE, subtitle="",
                      factor_levels=NULL) {

  # adapted from:
  # https://stackoverflow.com/questions/65269825/is-it-possible-to-recreate-the-functionality-of-bayesplots-mcmc-areas-plot-in
  # https://stackoverflow.com/questions/49961582/how-shade-area-under-ggridges-curve

  if(is.null(factor_levels)){
    factor_levels <- data %>% select(-n_sim) %>% names()
  }
  p_data <-
    data %>%
    pivot_longer(!n_sim)

  bw <- NULL

  p_data$name <- factor(p_data$name, levels=factor_levels)

  cred.int <- data %>%
    pivot_longer(!n_sim) %>%
    group_by(name) %>%
    summarise(CI=list(hdi(value, ci=ci)),
              m=median(value, na.rm=TRUE)) %>%
    unnest_wider(CI)

  p <- p_data %>%
    ggplot(aes(x = value, y = name)) +
    stat_density_ridges(
      geom = "density_ridges_gradient",
      calc_ecdf = TRUE,
      quantiles = c(0.5),
      scale=0.9,
      quantile_lines = TRUE,
      bandwidth = bw
    )

  d <- transform(ggplot_build(p)$data[[1]], name=group) %>%
    left_join(.,
              mutate(cred.int, name=as.numeric(factor(name, levels=factor_levels))),
              by="name") %>%
    filter(x < CI_high, x > CI_low)

  p <- p +
    geom_ribbon(
      data=d,
      aes(x, ymin=ymin, ymax=ymax, group=group),
      fill="#ADD8E6") +
    stat_density_ridges(
      geom = "density_ridges_gradient",
      calc_ecdf = TRUE,
      quantiles = c(0.5),
      scale=0.9,
      quantile_lines = TRUE,
      fill="transparent",
      bandwidth = bw
    ) +
    coord_flip()+
    theme_bw() +
    labs(
      y="", x="",
      subtitle=subtitle
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "position.none")

  if(limit_y){
    p <- p +
      scale_x_continuous(limits=c(0,1))
  } else {
    p <- p +
      scale_x_continuous(labels=scales::dollar_format())
  }
  p
}

plot_binned_ridges <- function(data, ci=0.95, hdi=T, limit_y=FALSE, subtitle="",
                               factor_levels=NULL) {

  p_data <-
    data %>%
    pivot_longer(!n_sim)

  if(hdi) {
    cred.int <- data %>%
      pivot_longer(!n_sim) %>%
      group_by(name) %>%
      summarise(CI=list(hdi(value, ci=ci)),
                m=median(value, na.rm=TRUE)) %>%
      unnest_wider(CI)

    p_data <-
      left_join(p_data, cred.int, by="name") %>%
      mutate(in_interval=value > CI_low & value < CI_high)
  } else {
    probs <- c((1 - ci)/2, 1 - (1 - ci)/2)

    p_data <-
      p_data %>%
      group_by(name) %>%
      arrange(value) %>%
      mutate(percentile=row_number()/n()) %>%
      ungroup() %>%
      mutate(in_interval = percentile > probs[1] & percentile < probs[2])
    # return(p_data)
  }


  if(is.null(factor_levels)){
    factor_levels <- data %>% select(-n_sim) %>% names()
  }

  p_data$name <- factor(p_data$name, levels=factor_levels)

  p <-
    p_data %>%
    ggplot(aes(x = value, y = name, height = stat(count), fill=in_interval)) +
    geom_density_ridges(stat = "binline", scale = 0.95, draw_baseline=FALSE) +
    coord_flip() +
    theme_bw() +
    labs(
      y="", x="",
      subtitle=subtitle
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "position.none") +
    scale_fill_manual(values=c("#ADD8E6","grey50"), breaks=c(TRUE, FALSE))

  if(!limit_y){
    p <- p + scale_x_continuous(labels=scales::dollar_format())
  }
  p
}

# extract content from lists (made from parallel processing of simulations)
extract_result_plots <- function(l) {
  res <- list()
  for(i in 1:length(l)) {
    res <- c(res, list(l[[i]]$plot_result))
  }
  res
}

extract_threshold_plots <- function(l) {
  res <- list()
  for(i in 1:length(l)) {
    res <- c(res, list(l[[i]]$plot_threshold))
  }
  res
}

extract_summaries <- function(l, format_max_nmb=T) {
  res <- list()
  for(i in 1:length(l)) {
    df_summary <- l[[i]]$summary
    if(format_max_nmb){
      medians <- as.numeric(str_extract(df_summary$summary, "-?\\d+\\.?\\d*"))
      df_summary$summary[which.max(medians)] <- paste0("<b>", df_summary$summary[which.max(medians)], "</br>")
    }
    res <- c(res, list(df_summary))
  }
  res
}

extract_thresholds_summary <- function(l) {
  res <- list()
  for(i in 1:length(l)) {
    res <- c(res, list(l[[i]]$thresholds_summary))
  }
  res
}
