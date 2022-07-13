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

summarize_sims <- function(x, prob, hdi, agg_fx, interval_in_text=F) {
  r <- function(x) round(x, digits=2)

  if(hdi){
    hdi <- hdi(x, ci=prob)
    res <- glue::glue(
      "{r(agg_fx(x))} [{r(hdi$CI_low)}, {r(hdi$CI_high)}]"
    )
  } else {
    probs <- c((1 - prob)/2, 1 - (1 - prob)/2)
    quantiles <- quantile(x, probs=probs)

    if(interval_in_text){
      res <- glue::glue(
        "{r(agg_fx(x))} [{scales::percent(prob)} Interval:{r(quantiles[[1]])}, {r(quantiles[[2]])}]"
      )
    }else {
      res <- glue::glue(
        "{r(agg_fx(x))} [{r(quantiles[[1]])}, {r(quantiles[[2]])}]"
      )
    }

  }
  list(summary=res)
}

get_summary <- function(data, sample_size, n_sims, n_valid, sim_auc, event_rate, inb_ref_col=NA,
                        agg_fx=median, hdi=F, ci=0.95, make_max_bold=T, recode_methods_vector=NULL, ...) {

  if(!is.na(inb_ref_col)) {
    data <-
      data %>%
      mutate(across(!n_sim, function(x) x - !!rlang::sym(inb_ref_col))) %>%
      select(-all_of(inb_ref_col))
  }

  df_summary <- data %>%
    pivot_longer(!n_sim, names_to="method", values_to="nmb")

  df_n_best <- df_summary %>%
    group_by(n_sim) %>%
    arrange(desc(nmb)) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(method) %>%
    summarize(n_best=n())


  df_summary <- as.data.table(df_summary)[
    ,
    summarize_sims(x=nmb, prob=ci, hdi=hdi, agg_fx=agg_fx),
    by = list(method)
  ]

  df_summary$sample_size <- sample_size
  df_summary$n_sims <- n_sims
  df_summary$n_valid <- n_valid
  df_summary$sim_auc <- sim_auc
  df_summary$event_rate <- event_rate

  df_summary <- left_join(df_summary, df_n_best, by="method") %>%
    mutate(n_best = ifelse(is.na(n_best), 0, n_best),
           percent_best = scales::percent(n_best/n_sims, accuracy=1),
           n_best_percent = glue::glue("{n_best} ({percent_best})"))

  if(make_max_bold){
    values <- as.numeric(str_extract(df_summary$summary, "-?\\d+\\.?\\d*"))
    df_summary$summary[values==max(values)] <- paste0("<b>", df_summary$summary[values==max(values)], "</b>")

    n_bests <- df_summary$n_best
    df_summary$n_best_percent[n_bests==max(n_bests)] <- paste0("<b>", df_summary$n_best_percent[n_bests==max(n_bests)], "</b>")
  }

  if(!is.null(recode_methods_vector)){
    recode_methods_vector
    df_summary$method <- recode_methods(df_summary$method, recode_methods_vector)
  }
  df_summary
}

recode_methods <- function(method, named_vector){
  recode(
    method,
    !!!setNames(names(named_vector), named_vector),
    .default=method
  )
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


get_plot_data <- function(data, factor_levels) {

  pivoted_data <- pivot_longer(data, !n_sim)

  if(is.null(factor_levels)){
    factor_levels <- data %>% select(-n_sim) %>% names()
  }

  pivoted_data$name <- factor(pivoted_data$name, levels=factor_levels)

  pivoted_data
}

add_interval <- function(data, ci=0.95, hdi=F) {
  if(hdi) {
    cred.int <-
      data %>%
      group_by(name) %>%
      summarise(CI=list(hdi(value, ci=ci)),
                m=median(value, na.rm=TRUE)) %>%
      unnest_wider(CI)

    data <-
      left_join(data, cred.int, by="name") %>%
      mutate(in_interval=value > CI_low & value < CI_high)
  } else {
    probs <- c((1 - ci)/2, 1 - (1 - ci)/2)

    data <-
      data %>%
      group_by(name) %>%
      arrange(value) %>%
      mutate(percentile=row_number()/n()) %>%
      ungroup() %>%
      mutate(in_interval = percentile > probs[1] & percentile < probs[2])
  }
  data
}

plot_binned_ridges <- function(data, ci=0.95, hdi=F, limit_y=FALSE, subtitle="",
                               factor_levels=NULL) {

  p_data <-
    get_plot_data(data, factor_levels=factor_levels) %>%
    add_interval(ci=ci, hdi=hdi)

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


plot_fw_histogram <- function(data, inb_ref_col=NA, ci=0.95, hdi=F, limit_y=FALSE, subtitle="",
                              factor_levels=NULL, agg_fx=median, n_bins=40,
                              n_breaks=3, plot_labels=labs(x="", y=""),
                              agg_line_alpha=0.6, agg_line_size=2, remove_axis=F,
                              label_wrap_width=10,
                              extra_theme=theme(panel.spacing  = unit(0, "lines"),
                                                axis.ticks.x = element_blank(),
                                                axis.text.x = element_blank(),
                                                strip.background = element_rect(fill="#f1f1f1")),
                              only_show_interval=F) {

  if(!is.na(inb_ref_col)) {
    data <-
      data %>%
      mutate(across(!n_sim, function(x) x - !!rlang::sym(inb_ref_col))) %>%
      select(-all_of(inb_ref_col))
  }

  p_data <-
    get_plot_data(data, factor_levels=factor_levels) %>%
    add_interval(ci=ci, hdi=hdi)

  df_agg <-
    p_data %>%
    group_by(name) %>%
    summarize(m=agg_fx(value))

  if(only_show_interval){
    p_data <- filter(p_data, in_interval)
    fill_cols <- c("grey50", "#ADD8E6")
  }else {
    fill_cols <- c("grey50", "grey50", "#ADD8E6")
  }

  p <-
    p_data %>%
    ggplot(aes(value, fill=in_interval)) +
    geom_histogram(bins=n_bins) +
    coord_flip() +
    facet_wrap(~name, labeller=label_wrap_gen(width=label_wrap_width), nrow=1) +
    theme_bw() +
    scale_fill_manual(values=c(fill_cols)) +
    guides(fill="none") +
    scale_y_continuous(n.breaks=n_breaks) +
    plot_labels

  if(!limit_y){
    p <- p + scale_x_continuous(labels=scales::dollar_format())
  }

  my_plot_innards <- ggplot_build(p)

  extracted_points <- tibble(
    outcome = my_plot_innards[["data"]][[1]][["x"]],
    count = my_plot_innards[["data"]][[1]][["y"]],
    in_interval = (my_plot_innards[["data"]][[1]][["group"]]) %>% as.factor(),
    method = (my_plot_innards[["data"]][[1]][["PANEL"]] %>% as.character)
  )

  heights <-
    df_agg %>%
    rownames_to_column(var="method_n") %>%
    left_join(extracted_points, by=c("method_n"="method")) %>%
    mutate(diff=abs(m-outcome)) %>%
    group_by(name) %>%
    arrange(diff) %>%
    slice(1)

  p <- p + geom_segment(data=heights, aes(x=m, xend=m, y=0, yend=count), size=agg_line_size, alpha=agg_line_alpha)

  if(remove_axis){
    p <- p + scale_y_continuous(breaks = NULL)
  }

  if(!is.null(extra_theme)){
    p <- p + extra_theme
  }

  p
}


# extract content from lists (made from parallel processing of simulations)
get_plot_list <- function(out_list,
                          rename_vector,
                          inb_ref_col=NA,
                          get_what = c("nmb", "inb", "cutpoints", "calibration"),
                          ...){

  get_what <- get_what[1]
  plotlist <- list()

  if(get_what=="calibration") {
    for(i in 1:length(out_list)){
      p <- out_list[[i]]$calibration_plot +
        scale_x_continuous(limits=c(0, 0.5)) +
        scale_y_continuous(limits=c(0, 0.5))
      plotlist <- c(plotlist, list(p))
    }
    return(plotlist)
  }



  for(i in 1:length(out_list)){
    if(get_what %in% c("nmb", "inb")){
      results <- out_list[[i]]$df_result
    } else {
      results <- out_list[[i]]$df_thresholds %>%
        select(where(~n_distinct(.) > 1))
    }
    if(!missing(rename_vector)){
      results <- rename(results, any_of(rename_vector))
    }
    p <- plot_fw_histogram(data=results, inb_ref_col=inb_ref_col, ...)
    plotlist <- c(plotlist, list(p))
  }
  plotlist
}


extract_summaries<- function(out_list,
                             rename_vector,
                             get_what = c("nmb", "inb", "cutpoints"),
                             inb_ref_col=NA,
                             agg_fx=median,
                             hdi=F,
                             ci=0.95,
                             ...) {

  get_what <- get_what[1]
  summarylist <- list()

  for(i in 1:length(out_list)){
    if(get_what %in% c("nmb", "inb")){
      results <- out_list[[i]]$df_result
    } else {
      results <- out_list[[i]]$df_thresholds
      results <- select(results, -treat_all, -treat_none)
    }
    if(!missing(rename_vector)){
      results <- rename(results, any_of(rename_vector))
    }
    summary_i <- do.call(
      get_summary,
      c(
        out_list[[i]]$meta_data,
        list(
          data=results,
          agg_fx=agg_fx,
          hdi=hdi,
          ci=ci,
          recode_methods_vector=rename_vector,
          inb_ref_col=inb_ref_col,
          make_max_bold=(get_what!="cutpoints")
        )
      )
    )
    summarylist <- c(summarylist, list(summary_i))
  }
  summarylist
}


make_table <- function(l, get_what = c("nmb", "inb", "cutpoints"),
                       rename_vector=cols_rename,
                       inb_ref_col=NA,
                       agg_fx=median, hdi=F, ci=0.95,
                       save_path=NULL) {
  tbl <- rbindlist(extract_summaries(
    l, rename_vector=rename_vector, get_what=get_what, agg_fx=agg_fx, hdi=hdi, ci=ci, inb_ref_col=inb_ref_col
  ))

  tbl <-
    tbl %>%
    group_by(sample_size, n_sims, n_valid, sim_auc, event_rate)  %>%
    mutate(.group_id=cur_group_id()) %>%
    ungroup() %>%
    select(-sample_size, -n_sims, -n_valid, -.group_id, -n_best, -percent_best, -n_best_percent) %>%
    select(Rate=event_rate, `Model AUC`=sim_auc, everything()) %>%
    pivot_wider(names_from="method", values_from="summary") %>%
    formattable() %>%
    kable(escape=F,
          caption=glue::glue("Data presented as median [{percent(ci, digits=0)} Intervals]")) %>%
    kable_styling()
  if(!is.null(save_path)) {
    save_kable(tbl, save_path)
  }
  tbl
}

keep_only_first_plot_strip <- function(plotlist){
  if(length(plotlist)==1) return(plotlist)
  for(i in 2:length(plotlist)){
    plotlist[[i]] <-
      plotlist[[i]] +
      theme(
        strip.text.x = element_blank()
      )
  }
  plotlist
}
