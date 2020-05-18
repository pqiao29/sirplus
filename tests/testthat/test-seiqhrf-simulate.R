test_that("seiqhrf and simulate_seiqhrf produce identical output", {
    
    s.num = 500
    q.num = 10
    nsteps = 10
    nsims = 3
    arec.rate = 0    ### has to be fixed 0 for comparison
    rec.rate = 0     ### has to be fixed 0 for comparison
    rec.rand = TRUE  ### has to be fixed TRUE for comparison
    
    No_seeds <- 10
    seed_list <- sample(1:1000, No_seeds)
    comp <- rep(NA, No_seeds)
    i <- 1
    for(seed in seed_list){
        Churhes_res <- simulate_seiqhrf(nsteps = nsteps, nsims = nsims, rec.rate = rec.rate, rec.rand = rec.rand, 
                                        arec.rate = arec.rate, s.num = s.num, q.num = q.num,
                                        infection.FUN = infection.seiqhrf.icm, 
                                        recovery.FUN = progress.seiqhrf.icm, 
                                        departures.FUN = departures.seiqhrf.icm, 
                                        arrivals.FUN = arrivals.seiqhrf.icm, seed = seed)$sim
        class(Churhes_res) <- "seiqhrf"
        
        param <- param_seiqhrf(arec.rate = arec.rate, rec.rate = rec.rate)
        init <- init_seiqhrf(s.num = s.num, q.num = q.num)
        control <- control_seiqhrf(nsteps = nsteps, nsims = nsims, rec.rand = rec.rand)
        sirplus_res <- seiqhrf(init, control, param, NULL, seed)
        
        comp[i] <- identical(Churhes_res[3:4], sirplus_res[3:4]) # Due to $usr.specified in control and param, can only compare "epi" and "times"
        i <- i + 1
    }
    
    expect_equal(sum(comp), No_seeds)
})
