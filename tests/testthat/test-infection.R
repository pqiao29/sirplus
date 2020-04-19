test_that("Identical output as Churches' original function: infection.FUN", {
    
    control <- control_seiqhrf()
    param <- param_seiqhrf()
    init <- init_seiqhrf()
    
    dat <- do.call(sirplus::initialize.FUN,  list(param, init, control))
    sd_No <- 10   ## Number of seeds to be tested
    comp <- rep(NA, sd_No)
    for(seed in 1:sd_No){
        dat1 <- do.call(infection.seiqhrf.icm, list(dat, at = 2, seed)) 
        dat2 <- do.call(infection.FUN, list(dat, at = 2, seed))         
        
        comp[seed] <- identical(dat1, dat2)
       
    }
    
    expect_equal(sum(comp), sd_No)
})