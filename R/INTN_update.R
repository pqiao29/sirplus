# Needed for recovery.FUN --------------------------------------------------------------
update_status <- function(rate, rand, active, status, label, state, at, expTime = NULL, prog.dist.scale = NULL, prog.dist.shape = NULL){
    
    smp_sz <- 0
    at_idx <- NULL
    
    idsElig <- which(active == 1 & status %in% label)
    nElig <- length(idsElig)
    
    if (nElig > 0) {

        if (rand) {
            vecProg <- which(stats::rbinom(nElig, 1, rate) == 1)
            if (length(vecProg) > 0) {
                at_idx <- idsElig[vecProg]
                smp_sz <- length(at_idx)
                status[at_idx] <- state
            }
        }else{
            vecTimeSinceExp <- at - expTime[idsElig]
            vecTimeSinceExp[is.na(vecTimeSinceExp)] <- 0
            gammaRatesElig <- stats::pweibull(vecTimeSinceExp, prog.dist.shape, scale=prog.dist.scale)
            smp_sz <- round(sum(gammaRatesElig, na.rm=TRUE)) 
            smp_prob <- gammaRatesElig
            
            if(smp_sz > 0){
                at_idx <- EpiModel::ssample(idsElig, smp_sz, prob = smp_prob)
                status[at_idx] <- state
            }
            
        }
        
    }
    
    list(smp_sz, at_idx, status)
}

# Needed for departures.FUN ------------------------------------------------------------
update_active <- function(rate, rand, active, status, label){
    
    smp_sz <- 0
    act_idx <- NULL
    
    idsElig <- which(active == 1 & status == label)
    nElig <- length(idsElig)
    
    if (nElig > 0) {
        
        gElig <- rep(1, nElig)
        
        if (rand) {
            vec <- which(stats::rbinom(nElig, 1, rate) == 1)
            
            if (length(vec) > 0) {
                act_idx <- idsElig[vec]
                smp_sz <- length(act_idx)
            }
        } else {
            smp_sz <- min(round(sum(rate[gElig == 1])), sum(gElig == 1))
            act_idx <- EpiModel::ssample(idsElig[gElig == 1], smp_sz)
        }
    }
    
    return(list(smp_sz, act_idx))
}