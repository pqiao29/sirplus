test_that("Identical output as Churches' original function: recovery.FUN", {
    
    control <- control_seiqhrf(rec.rand = TRUE)
    param <- param_seiqhrf(arec.rate = 0, rec.rate = 0)
    init <- init_seiqhrf()
    
    at <- 2
    dat <- do.call(initialize.FUN, list(param, init, control))
    dat <- do.call(infection.FUN, list(dat, at))   
    
    No_seeds <- 10
    seed_list <- sample(1:1000, No_seeds)
    comp <- rep(NA, No_seeds)
    i <- 1
    for(seed in seed_list){
        dat1 <- do.call(recovery.FUN, list(dat, at, seed))
        dat2 <- do.call(progress.seiqhrf.icm, list(dat, at, seed))
        comp[i] <- identical(dat1, dat2)
        i <- i + 1
    }
    
    expect_equal(sum(comp), No_seeds)
})