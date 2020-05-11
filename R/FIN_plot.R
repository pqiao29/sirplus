#' Plot simulation result 
#'
#' @param x An seiqhrf object returned from function \code{\link{seiqhrf}}.
#' @param method If "times", plot Duration frequency distributions.
#'               If "weekly_local", plot local weekly estimates from simulation.
#'               If NULL, plot sirplus plots.
#' @param return_df In effect only when method == "weekly", if TRUE returns 
#'        also the dataframe used for plotting as well as the ggplot object.
#' @param comp_remove Compartments to remove. Suggest c(s.num, r.num)
#' @param time_limit Number of steps (days) to plot. 
#' @param ci T/F to include 95\% confidence intervals in sirplus plot.
#' @param sep_compartments T/F use faceting to show each compartment in a 
#'        separate plot, only works if plotting a single simulation.
#' @param trans Y-axis transformation (e.g. log2, log10). 
#' @param known Dataframe with known compartment numbers to plot alongside
#'        projections
#' @param start_date Date for day 0. 
#' @param show_start_date First date to show in plots. Use ymd format. If FALSE,
#'        shows from step 1. 
#' @param x_axis Title for x-axis. 
#' @param plot_title Title for whole plot. 
#' @param market.share between 0 and 1, percentage of local hospital beds in 
#'        the simulated unit (e.g. state)
#' @param icu_percent between 0 and 1, percentage of patients that should go to 
#'        ICU among the ones that need hospitalization
#' @param sim_population Size of population simulated. Only needed if providing 
#'        `total_population`.
#' @param total_population True population size, needed only if simulation size 
#'        is smaller than the true population size due to computational cost 
#'        etc.
#' @param ... Additional parameters
#' 
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom dplyr "%>%"
#' @importFrom dplyr bind_rows
#' @importFrom dplyr filter
#' @importFrom lubridate ymd
#' @import ggplot2
#' @export
plot.seiqhrf <- function(x, method = NULL, 
                         comp_remove = "none",
                         time_limit = 90,
                         ci = TRUE,
                         sep_compartments = FALSE,
                         trans = FALSE,
                         known = NULL,
                         start_date = ymd("2020-03-21"),
                         show_start_date = FALSE,
                         x_axis = 'Date (MM-DD)',
                         plot_title = '', 
                         return_df = TRUE, 
                         market.share = .04,
                         icu_percent = .1, 
                         sim_population = 1000,
                         total_population = NULL, ...) {
    
    if(is.null(method)){
        plot_sirplus(x, comp_remove = comp_remove,
                     time_limit = time_limit,
                     ci = ci,
                     sep_compartments = sep_compartments,
                     trans = trans,
                     known = known,
                     start_date = start_date,
                     x_axis = x_axis,
                     plot_title = plot_title,
                     sim_population = sim_population,
                     total_population = total_population)
        
    }else if(method == "times"){
        
            if(!inherits(x, "seiqhrf")) stop(
                "If method == times, x needs to be an seiqhrf object")
            plot_times(x) 
            
    }else if(method == "weekly_local"){
        
            ret <- get_weekly_local(x, market.share = market.share,
                                    icu_percent = icu_percent, 
                                    start_date = start_date,
                                    show_start_date = show_start_date,
                                    time_limit = time_limit,
                                    sim_population = sim_population,
                                    total_population = total_population)
             if(return_df){
                 return(ret) 
             }else{
                 return(ret$plot)
             }
    }
}

#' Plot simulation result given list
#'
#' @param x A list of seiqhrf objects returned from \code{\link{seiqhrf}}.
#' @param comp_remove Compartments to remove. Suggest c(s.num, r.num)
#' @param time_limit Number of steps (days) to plot. 
#' @param ci T/F to include 95\% confidence intervals in sirplus plot.
#' @param sep_compartments T/F use faceting to show each compartment in a 
#'        separate plot, only works if plotting a single simulation.
#' @param trans Y-axis transformation (e.g. log2, log10).  
#' @param known Dataframe with known compartment numbers to plot alongside
#'        projections
#' @param start_date Date for day 0. 
#' @param x_axis Title for x-axis. 
#' @param plot_title Title for whole plot. 
#' total_population
#' @param sim_population Size of population simulated. Only needed if providing 
#'        `total_population`.
#' @param total_population True population size, needed only if simulation size 
#'        is smaller than the true population size due to computational cost 
#'        etc.
#' @param ... Additional parameters
#' 
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom dplyr "%>%"
#' @importFrom dplyr bind_rows
#' @importFrom dplyr filter
#' @importFrom lubridate ymd
#' @import ggplot2
#' @export
plot.list <- function(x, comp_remove = "none",
                      time_limit = 90,
                      ci = TRUE,
                      sep_compartments = FALSE,
                      trans = FALSE,
                      known = NULL,
                      start_date = ymd("2020-03-21"),
                      x_axis = 'Date (MM-DD)',
                      plot_title = '',
                      sim_population = 1000,
                      total_population = NULL, ...){
    
    plot_sirplus(x, comp_remove = comp_remove,
                 time_limit = time_limit,
                 ci = ci,
                 sep_compartments = sep_compartments,
                 trans = trans,
                 known = known,
                 start_date = start_date,
                 x_axis = x_axis,
                 plot_title = plot_title,
                 sim_population = sim_population,
                 total_population = total_population)
}




#' Wrapper for primary sirplus plotting function
#'
#' Flexible function to generate sirplus plots (i.e. compartment counts over 
#' time). This function allows for plotting multiple experiments, viewing the 
#' plots of different scales (e.g. log2), plotting compartments separately,
#' adding 95\% CIs, and plotting known data along side the simulations.  
#'
#' @param x A seiqhrf object (or list of multiple seiqhrf objects) returned 
#'        from \code{\link{seiqhrf}}.
#' @param comp_remove Compartments to remove. Suggest c(s.num, r.num)
#' @param time_limit Number of steps (days) to plot. 
#' @param ci T/F to include 95\% confidence intervals in sirplus plot.
#' @param sep_compartments T/F use faceting to show each compartment in a 
#'        separate plot, only works if plotting a single simulation.
#' @param trans Y-axis transformation (e.g. log2, log10). 
#' @param known Dataframe with known compartment numbers to plot alongside
#'        projections
#' @param start_date Date for day 0. 
#' @param x_axis Title for x-axis. 
#' @param plot_title Title for whole plot. 
#' @param sim_population Size of population simulated. Only needed if providing 
#'        `total_population`.
#' @param total_population True population size, needed only if simulation size 
#'        is smaller than the true population size due to computational cost 
#' @param ... Additional parameters
#' 
#' @return ggplot2 object
#' 
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom dplyr "%>%"
#' @importFrom dplyr bind_rows
#' @importFrom dplyr filter
#' @import ggplot2
#' @export
plot_sirplus <- function(x, comp_remove,
                         time_limit,
                         ci,
                         sep_compartments,
                         trans,
                         known,
                         start_date,
                         x_axis,
                         plot_title,
                         sim_population,
                         total_population,...){
    
    # Convert from seiqhrf object to dataframe
    plot_df <- format_sims(x, time_limit = time_limit, start_date = start_date)
    reo_exp <- function(x) {factor(x, levels = unique(plot_df$experiment))}
    
    # Get Confidence Intervals
    if(ci){
        plot_df <- get_ci(x, plot_df)
    }
    
    # Scale up to full population size if needed
    if(!is.null(total_population)){
      if(total_population < sim_population) 
        stop("total population should be larger than simulated size")
      
      scale_factor <- total_population/sim_population
      plot_df <- plot_df %>% mutate(count = count * scale_factor)
      
      if(ci){
        plot_df <- plot_df %>% mutate(CI.1 = count * scale_factor,
                                      CI.2 = CI.2 * scale_factor,
                                      qntCI.1 = qntCI.1 * scale_factor,
                                      qntCI.2 = qntCI.2 * scale_factor)
      }
    } 
    
    # Add known compartment counts
    if(is.data.frame(known)){
        plot_df <- add_known(plot_df, known = known, start_date = start_date)
    }
    
    # Define compartment names and colours
    comps <- c("s.num", "e.num", "i.num", "q.num", "h.num", "r.num", "f.num")
    compcols <- c(s.num = "#4477AA", e.num = "#66CCEE", i.num = "#CCBB44", 
                  q.num = "#AA3377", h.num = "#EE6677", r.num = "#228833", 
                  f.num = "#BBBBBB")
    complabels <- c(s.num = "S: Susceptible", e.num = "E: Asymptomatic", 
                    i.num = "I: Infected", q.num = "Q: Self-isolated", 
                    h.num = "H: Hospitalized", r.num = "R: Recovered",
                    f.num = "F: Case Fatalities")
    
    # Filter compartments
    comp_plot <- setdiff(comps, comp_remove)
    plot_df <- plot_df %>% filter(compartment %in% c(comp_plot)) %>%
      filter(time <= time_limit)
  
    # Plot with options
    p <- ggplot(plot_df, aes(x = Date, y = count, colour = compartment, 
                             linetype = sim)) + 
        geom_line(size = 1.2, alpha = 0.8) + 
        scale_x_date(date_breaks = "1 week", date_labels = "%m-%d") + 
        scale_colour_manual(values = compcols, labels = complabels) + 
        labs(title = plot_title, x = x_axis, y = "Prevalence") +
        theme_bw() + theme(axis.text.x = element_text(angle = 90))
    
    if(length(unique(plot_df$experiment)) > 1){
        p <- p + facet_grid(reo_exp(experiment) ~ ., scales = 'free')
    }
    
    if(sep_compartments){
        p <- p + facet_grid(compartment ~ ., scales = 'free')
    }
    
    if(trans){
        p <- p + scale_y_continuous(trans = trans) 
    }
    
    if(ci){
        p <- p + geom_ribbon(aes(ymin=qntCI.1, ymax=qntCI.2, x=Date, 
                                 fill = compartment, colour = NULL),
                             alpha = 0.4) +
            scale_color_manual(values = compcols, labels = complabels) +
            scale_fill_manual(values = compcols, guide = FALSE)
    }
    
    p
}

#' Plot compartment duration distributions
#'
#' Function to plot Duration frequency distributions. If multiple simulations 
#' were performed (nsim >1), durations from sims are appended to each other.
#'
#' @param sim An seiqhrf object returned from function \code{\link{seiqhrf}}.
#' 
#' @return ggplot2 object
#' 
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom dplyr "%>%"
#' @importFrom dplyr bind_rows
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' 
#' @export
plot_times <- function(sim) {
    
    for (s in 1:sim$control$nsims) {
        if (s == 1) {
            times <- sim$times[[paste("sim", s, sep = "")]]
            times <- times %>% mutate(s = s)
        } else {
            times <- times %>% bind_rows(sim$times[[paste("sim", s, sep = "")]] 
                                         %>% mutate(s = s))
        }
    }
    
    times <- times %>% mutate(infTime = ifelse(infTime < 0, -5, infTime), 
                              expTime = ifelse(expTime < 0, -5, expTime)) %>% 
        mutate(incubation_period = infTime - expTime, 
               illness_duration = recovTime - expTime, 
               illness_duration_hosp = dischTime - expTime, 
               hosp_los = dischTime - hospTime, 
               quarantine_delay = quarTime - infTime,
               survival_time = fatTime - infTime) %>% 
        select(s, incubation_period, quarantine_delay, illness_duration, 
               illness_duration_hosp, hosp_los, survival_time) %>% 
        pivot_longer(-s, names_to = "period_type", values_to = "duration") %>% 
        mutate(period_type = factor(period_type, 
                                    levels = c("incubation_period", 
                                               "quarantine_delay", 
                                               "illness_duration", 
                                               "illness_duration_hosp", 
                                               "hosp_los", 
                                               "survival_time"), 
                                    labels = c("Incubation\nperiod", 
                                               "Delay entering\nisolation", 
                                               "Illness\nduration",
                                               "Illness\nduration (hosp)", 
                                               "Hospital stay\nduration", 
                                               "Survival time\nfor fatalities"),
                                    ordered = TRUE))
    
    times %>% filter(duration <= 30) %>% ggplot(aes(x = duration)) + 
        geom_bar() + facet_grid(period_type ~ ., scales = "free_y") + 
        labs(title = "Compartment Duration Distributions")
}


#' Extract and plot information of local and weekly estimates from simulation 
#'
#' @param sim An \code{seiqhrf} object returned by \link{simulate_seiqhrf}. 
#' @param market.share between 0 and 1, percentage of local hospital beds in 
#'        the simulated unit (e.g. state)
#' @param icu_percent between 0 and 1, percentage of patients that should go to 
#'        ICU among the ones that need hospitalization
#' @param start_date Epidemic start date. Default is 'na', if not provided will 
#'        plot week numbers, if provided will plot the first day (Sunday) of the
#'        week.
#' @param show_start_date First date to show in plots. Use ymd format. If FALSE,
#'        shows from step 1.
#' @param time_limit Number of days to include. 
#' @param sim_population Size of population simulated. Only needed if providing 
#'        `total_population`.
#' @param total_population True population size, needed only if simulation size 
#'        is smaller than the true population size due to computational cost 
#'        etc.
#' 
#' @return 
#' \itemize{
#' \item \code{plot:} A \code{ggplot} object, bar charts of count of patients 
#'              requiring hospitalization and ICU respectively
#' \item \code{result:} A dataframe
#'    \itemize{\item \code{week:}  week number from input \code{sim},
#'             \item \code{hosp:} the number of patients that require hospitalization locally,
#'             \item \code{icu:} the number of patients that require ICU locally. }
#
#' }
#' 
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom utils head
#' @export
get_weekly_local <- function(sim, 
                             market.share,
                             icu_percent, 
                             start_date,
                             show_start_date,
                             time_limit,
                             sim_population,
                             total_population){

    # Get h.num and 95% quantile CIs
    sim_mean <- as.data.frame(sim, out = "mean")
    ci_info <- as.data.frame.list(summary.seiqhrf(sim)$h.num)
    hosp <- data.frame('h.num' = sim_mean$h.num, 'ci5' = ci_info$qntCI.1,
                       'ci95' = ci_info$qntCI.2)
    hosp[is.na(hosp)] <- 0
    hosp <- hosp[1: time_limit, ]
    hosp$date <- start_date + as.numeric(row.names(hosp))
    
    if(class(show_start_date) == "Date"){
      hosp <- hosp %>% filter(date >= show_start_date)
    }
    
    # Scale for population size and hospital market share if needed
    if(!is.null(total_population)){
        if(total_population < sim_population) 
            stop("total Population should be larger than simulated size")
        cat("Scalling w.r.t total population")
        date_tmp <- hosp$date
        hosp$date <- NULL
        print(utils::head(hosp))
        hosp <- hosp*total_population/sim_population
        hosp$date <- date_tmp
    } 
    
    if(market.share < 0 || market.share > 1) stop("Market share has to be between 
                                                0 and 1")
    if(icu_percent < 0 || icu_percent > 1) stop("ICU percentage has to be between
                                              0 and 1")
    
    # Get weekly sums & calculate projected icu numbers
    hosp.wk <- hosp %>% group_by(yr_wk = cut(date, "week", start.on.monday = TRUE)) %>% 
        summarise(h.num = sum(h.num), h.ci5 = sum(ci5), h.ci95=sum(ci95)) %>%
        mutate(icu.num = h.num * icu_percent,
               icu.ci5 = h.ci5 * icu_percent,
               icu.ci95 = h.ci95 * icu_percent) %>%
        mutate(h.num = h.num - icu.num,
               h.ci5 = h.ci5 - icu.ci5,
               h.ci95 = h.ci95 - icu.ci95) 
    
    # Make long format for ggplot
    hosp.wk2 <- hosp.wk %>% pivot_longer(-yr_wk, names_to = c('type', 'metric'), 
                     values_to = 'val',
                     names_pattern = '(h|icu).(num|ci5|ci95)') %>%
        pivot_wider(names_from = metric, values_from = val)

    p <- ggplot(hosp.wk2, aes(x = yr_wk, y = num, color = type)) +
        geom_point() + 
        labs(y="Weekly Cumulative Count", x = "Week Start (Monday: YYYY-MM-DD)") +
        geom_errorbar(aes(ymin=ci5, ymax=ci95), size=0.5, width=.25) +
        scale_color_discrete(name = "Type", labels = c("General", "ICU")) + 
        theme_bw() + theme(axis.text.x = element_text(angle = 90))
    
    return(list("plot" = p, "result" = hosp.wk))
    
}



#' Format seiqhrf objects into dataframe for ggplot
#'
#' @param x An seiqhrf object returned from function \code{\link{seiqhrf}}.
#' @param time_limit Number of steps (days) to plot. 
#' @param start_date Date for day 0. 
#' 
#' @return dataframe
#' 
#' @export
format_sims <- function(x, time_limit, start_date){
    
    # Merge models to plot together
    if(class(x) == "seiqhrf"){
        sim_id <- "seiqhrf model" 
        plot_df <- as.data.frame(x, out = "mean")
        plot_df <- plot_df %>% mutate(experiment = sim_id)
    }else{
        
        sim_id <- names(x)
        if(is.null(sim_id)) stop("Please assign a name to each element in sims")
        
        plot_df <- as.data.frame(x[[1]], out = "mean")
        plot_df <- plot_df %>% mutate(experiment = sim_id[1])
        if(length(sim_id) > 1){
            for (i in (2:length(sim_id))) {
                tmp_df <- as.data.frame(x[[i]], out = "mean")
                tmp_df <- tmp_df %>% mutate(experiment = sim_id[i])
                plot_df <- plot_df %>% bind_rows(plot_df, tmp_df)
            }
        }
    }
    
    plot_df <- plot_df %>% filter(time <= time_limit) %>% 
        pivot_longer(-c(time, experiment),
                     names_to = "compartment", 
                     values_to = "count", 
                     values_ptypes = list(compartment = 'character',
                                          count = numeric())) %>%
        mutate(sim = 'sim',
               Date = start_date + time)

    return(plot_df)
}

#' Internal
#' @param obj An seiqhrf object
#' @param exp_name Name of experiment 

ci_info_update <- function(obj, exp_name = "seiqhrf model"){  
  
  ci_tmp <- as.data.frame.list(summary(obj))
  ci_tmp <- ci_tmp %>% 
    mutate(time = as.numeric(row.names(ci_tmp))) %>%
    pivot_longer(cols = -time, names_to = 'compartment',
                 values_to = 'mean',
                 values_ptypes = list(compartment = 'character',
                                      mean = numeric())) %>%
    tidyr::separate(compartment, 
                    into = c('compartment', 'metric'), 
                    sep='num.') %>%
    mutate(compartment = paste0(compartment, 'num'),
           experiment = exp_name) %>%
    pivot_wider(names_from = metric, values_from = mean) 
  
  ci_tmp
}

#' Get 95\% confidence intervals
#'
#' @param x An seiqhrf object returned from function \code{\link{seiqhrf}} or a list of seiqhrf objects.
#' @param plot_df Dataframe with known compartment numbers to plot alongside
#'        projections
#' 
#' @return dataframe with CIs and sd added
#' @importFrom tidyr separate
#' @importFrom dplyr full_join
#' 
#' @export
#' 
get_ci <- function(x, plot_df){
  
  # Get sim variance metrics for single seiqhrf object
  if(class(x) == "seiqhrf"){
    ci_info <- ci_info_update(x)
  }else{
    
    sim_id <- names(x)
    if(class(x) != "list") stop("The class of x should either be list or seiqhrf")
    ci_info <- ci_info_update(x[[1]], sim_id[1])
    
    if(length(sim_id) > 1){
      for (i in (2:length(sim_id))) {
        ci_tmp <- ci_info_update(x[[i]], sim_id[i])
        ci_info <- ci_info %>% bind_rows(ci_info, ci_tmp)
      }
    }
  }
  
  ci_info <- ci_info %>% mutate(sim = 'sim') %>%
    mutate(sim = 'sim')
  ci_info[is.na(ci_info)] <- 0
  
  plot_df <- plot_df %>% 
    dplyr::full_join(ci_info, by = c('time', 'compartment', 'experiment', 'sim'))
  
  return(plot_df)
}

#' Add known counts to sims dataframe for ggplot
#'
#' @param plot_df An seiqhrf object returned from function \code{\link{seiqhrf}}.
#' @param known Dataframe with known compartment numbers to plot alongside
#'        projections
#' @param start_date Date for day 0. 
#' 
#' @return dataframe with known data added
#' 
#' @export
#' 
add_known <- function(plot_df, known, start_date){
    
    # Add Date to known data
    missing_cols <- setdiff(names(plot_df), names(known))
    
    if("Date" %in% missing_cols){
        known$Date = start_date + known$time
    }
    
    known <- known %>% pivot_longer(cols = -c(time, Date), 
                                    names_to = 'compartment',
                                    values_to = 'count',
                                    values_ptypes = list(compartment = 'character',
                                                         count = numeric()))
    exps <- unique(plot_df$experiment)
    for(i in exps){
        known <- known %>% mutate(experiment = i, sim = 'known')
        missing_cols <- setdiff(names(plot_df), names(known))

        if(length(missing_cols) > 0){
            add <- rep(0, length(missing_cols))
            names(add) <- missing_cols
            known <- known %>% mutate(!!! add)
        }
        plot_df <- rbind(plot_df, known)
    } 
    
    return(plot_df)
}
  

