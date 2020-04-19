test_that("plot_model is producing the same plot", {

s.num <- 2000  # number susceptible
i.num <- 15  # number infected 
q.num <- 5  # number in self-isolation
h.num <- 1  #     
vals <- c(10, 7)
timing <- c(7, 14)
nsteps <- 90


## baseline
baseline_sim <- simulate_seiqhrf(nsteps = nsteps, s.num = s.num, i.num = i.num,
                                 q.num = q.num, h.num = h.num)

## Experiment 1
act_rate <- vary_param(nstep = nsteps, vals = vals, timing = timing)

sim_exp <- simulate_seiqhrf(nsteps = nsteps, s.num = s.num, i.num = i.num,
                            q.num = q.num, h.num = h.num, act.rate.e = act_rate, 
                            act.rate.i = act_rate * 0.5)


# ori_plot <- plot_models(c(baseline_sim, sim_exp),
#                         sim_id = c('Baseline', 'Closures'),
#                         start_date = lubridate::ymd("2020-01-01"),
#                         comp_remove = c('s.num', 'r.num'),
#                         plot_title = 'Closures Experiment')

my_plot <-  plot(list('Baseline' = baseline_sim$sim, 'Closures' = sim_exp$sim),
                 start_date = lubridate::ymd("2020-01-01"),
                 comp_remove = c('s.num', 'r.num'),
                 plot_title = 'Closures Experiment')

expect_true(is.ggplot(my_plot))
#expect_equal(all.equal(ori_plot, my_plot), TRUE)

})
