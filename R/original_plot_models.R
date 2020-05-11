#' Plot models
#'
#' Function to plot individuals models or multiple models for comparison.
#'
#' @param sims Single or list of sims to plot
#' @param sim_id String or list of strings to name each facet
#' @param comp_remove Compartments to remove. Suggest c(s.num, r.num)
#' @param time_lim Number of steps (days) to plot. 
#' @param trans Y-axis transformation (e.g. log2, log10). Default = none. 
#' @param known Dataframe with known compartment numbers to plot alongside
#'        projections
#' @param start_date Date for day 0. Default: ymd("2020-03-21"),
#' @param x_axis Title for x-axis. Default: 'Days since beginning of epidemic'
#' @param plot_title Title for whole plot. Default: 'SEIQHRF plot'
#'
#' @return ggplot2 object
#' 
#' @importFrom tidyr pivot_longer

plot_models <- function(sims,
                        sim_id = "baseline",
                        comp_remove = "none",
                        time_lim = 100,
                        trans = 'na',
                        known = NULL,
                        start_date = ymd("2020-03-21"),
                        x_axis = 'Days since beginning of epidemic',
                        plot_title = 'SEIQHRF'){
    
    # Define a standard set of colours to represent compartments
    comps <- c("s.num", "e.num", "i.num", "q.num", "h.num", "r.num", "f.num")
    compcols <- c(s.num = "#4477AA", e.num = "#66CCEE", i.num = "#CCBB44", 
                  q.num = "#AA3377", h.num = "#EE6677", r.num = "#228833", 
                  f.num = "#BBBBBB")
    complabels <- c(s.num = "S: Susceptible", e.num = "E: Asymptomatic", 
                    i.num = "I: Infected", q.num = "Q: Self-isolated", 
                    h.num = "H: Hospitalized", r.num = "R: Recovered",
                    f.num = "F: Case Fatalities")
    
    # Merge models to plot together
    for (i in (1: length(sim_id))) {
        if (i == 1) {
            plot_df <- sims[i*2]$df
            plot_df <- plot_df %>% mutate(experiment = sim_id[i])
        } else{
            tmp_df <- sims[i*2]$df
            tmp_df <- tmp_df %>% mutate(experiment = sim_id[i])
            plot_df <- plot_df %>% bind_rows(plot_df, tmp_df)
        }
    }
    
    plot_df <- plot_df %>% filter(time <= time_lim) %>% 
        pivot_longer(-c(time, experiment),
                     names_to = "compartment", 
                     values_to = "count")
    plot_df$sim <- 'sim'
    plot_df$Date <- start_date + plot_df$time
    plot_df$time <- NULL
    
    # Add known data
    if(is.data.frame(known)){
        exps <- unique(plot_df$experiment)
        for(i in exps){
            known$experiment <- i
            plot_df <- rbind(plot_df, known)
        } 
    }
    
    # Filter compartments
    comp_plot <- setdiff(comps, comp_remove)
    plot_df <- plot_df %>% filter(compartment %in% c(comp_plot))

    reo_exp <- function(x) { factor(x, levels = sim_id)}
    
    
    # Plot single model
    if(trans == "na"){
        if(length(sim_id) == 1){
            plot_df %>% ggplot(aes(x = Date, y = count+1, colour = compartment, 
                                   linetype = sim)) + 
                geom_line(size = 1.5, alpha = 0.8) + 
                scale_x_date(date_breaks = "1 week", date_labels = "%m-%d") + 
                scale_colour_manual(values = compcols, labels = complabels) + 
                labs(title = plot_title, x = x_axis, y = "Prevalence") +
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 90))
        } else (
            plot_df %>% ggplot(aes(x = Date, y = count, colour = compartment,
                                   linetype = sim)) + 
                facet_grid(reo_exp(experiment) ~ ., scales = 'free') + 
                scale_x_date(date_breaks = "1 week", date_labels = "%m-%d") + 
                geom_line(size = 1.5, alpha = 0.8) + 
                scale_colour_manual(values = compcols, labels = complabels) + 
                labs(title = plot_title, x = x_axis, y = "Prevalence") +
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 90))
        )
    }else{
        if(length(sim_id) == 1){
            plot_df %>% ggplot(aes(x = Date, y = count+1, colour = compartment, 
                                   linetype = sim)) + 
                scale_y_continuous(trans = trans) +
                scale_x_date(date_breaks = "1 week", date_labels = "%m-%d") + 
                geom_line(size = 1.5, alpha = 0.8) + 
                scale_colour_manual(values = compcols, labels = complabels) + 
                labs(title = plot_title, x = x_axis, y = "Prevalence") +
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 90))
        } else (
            plot_df %>% ggplot(aes(x = Date, y = count+1, colour = compartment, 
                                   linetype = sim)) + 
                facet_grid(reo_exp(experiment) ~ ., scales = 'free') + 
                scale_y_continuous(trans = trans) +
                scale_x_date(date_breaks = "1 week", date_labels = "%m-%d") + 
                geom_line(size = 1.5, alpha = 0.8) + 
                scale_colour_manual(values = compcols, labels = complabels) + 
                labs(title = plot_title, x = x_axis, y = "Prevalence") +
                theme_bw() + 
                theme(axis.text.x = element_text(angle = 90))
        )
    }
}
