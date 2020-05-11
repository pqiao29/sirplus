utils::globalVariables(c("time", "Date", "experiment", "compartment", "metric", 
                         "h.num", "ci5", "ci95", "h.ci5", "h.ci95", "icu.num", "icu.ci5", "icu.ci95", 
                         "yr_wk", "val", "num", "type",    # get_ci
                         "count", "sim",  # plot_models (original scripts)
                         "CI.1", "CI.2", "qntCI.1", "qntCI.2", # plot_sirplus
                         "infTime", "expTime", "recovTime", "dischTime","hospTime", "quarTime", "fatTime",
                         "incubation_period", "quarantine_delay", "illness_duration", "illness_duration_hosp", 
                         "hosp_los", "survival_time", "period_type", "duration" # plot_times
                         ))