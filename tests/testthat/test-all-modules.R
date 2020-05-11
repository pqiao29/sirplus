test_that("All modules(MOD) produce the same result as Churches' modules (.seiqhrf.icm)", {
    
    param <- param_seiqhrf(arec.rate = 0, rec.rate = 0)
    init <- init_seiqhrf(s.num = 1000)
    
    #### default functions: initialize.FUN, infection.FUN, recovery.FUN, departures.FUN, arrivals.FUN
    control1 <- control_seiqhrf(nsteps = 10, rec.rand = TRUE)
    ### Churches' original function:
    control2 <- control1
    control2$initialize.FUN <- "initialize.icm"
    control2$infection.FUN <- "infection.seiqhrf.icm"
    control2$recovery.FUN <- "progress.seiqhrf.icm"
    control2$departures.FUN <- "departures.seiqhrf.icm"
    control2$arrivals.FUN <- "arrivals.seiqhrf.icm"
    
    No_seeds <- 10
    seed_list <- sample(1:1000, No_seeds)
    comp <- rep(NA, No_seeds)
    i <- 1
    for(seed in seed_list){
        sim1 <- seiqhrf(param = param, init = init, control = control1, seed)
        sim2 <- seiqhrf(param = param, init = init, control = control2, seed)
        comp[i] <- identical(sim1, sim2)
        i <- i + 1
    }
    
    expect_equal(sum(comp), No_seeds)
})
