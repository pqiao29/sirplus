test_that("plot_model is producing the same plot", {

s.num <- 2000  # number susceptible
i.num <- 15  # number infected 
q.num <- 5  # number in self-isolation
h.num <- 1  #     
vals <- c(10, 7)
timing <- c(7, 14)
nsteps <- 90


control <- control_seiqhrf(nsteps = nsteps)
param <- param_seiqhrf()
init <- init_seiqhrf(s.num = s.num, i.num = i.num, q.num = q.num, h.num = h.num)

baseline_sim <- seiqhrf(init, control, param)

## Experiment 1
act_rate <- vary_param(nstep = nsteps, vals = vals, timing = timing)
param_exp <- param_seiqhrf(act.rate.e = act_rate, act.rate.i = act_rate * 0.5)
sim_exp <- seiqhrf(init, control, param)

exp_plot <-  plot(list('Baseline' = baseline_sim, 'Closures' = sim_exp),
                 start_date = lubridate::ymd("2020-01-01"),
                 comp_remove = c('s.num', 'r.num'),
                 plot_title = 'Closures Experiment')

expect_true(ggplot2::is.ggplot(exp_plot))

})
