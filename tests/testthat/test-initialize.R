test_that("Identical output as Churches' original function: initialize.FUN", {
    
    control <- control_seiqhrf()
    param <- param_seiqhrf()
    init <- init_seiqhrf()
    
    No_seeds <- 10
    seed_list <- sample(1:1000, No_seeds)
    comp <- rep(NA, No_seeds)
    i <- 1
    for(seed in seed_list){
        dat1 <- do.call(initialize.icm, list(param, init, control, seed))
        dat2 <- do.call(initialize.FUN, list(param, init, control, seed))
        comp[i] <- identical(dat1, dat2)
        i <- i + 1
    }
    
    expect_equal(sum(comp), No_seeds)
})
