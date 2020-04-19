## internal function

get_del <- function(dat, lab, acts, letgo = FALSE, seed = NULL){
    
    #if(!is.null(seed)) set.seed(seed)
    
    ## Edgelist
    p1 <- EpiModel::ssample(which(dat$attr$active == 1 & dat$attr$status != "f"), acts, replace = TRUE)  ## Bootstrap index, allow repetitive elements
    p2 <- EpiModel::ssample(which(dat$attr$active == 1 & dat$attr$status != "f"), acts, replace = TRUE)
    
    #### keep samplling until all elements in p1 and p2 differ
    del <- NULL
    if (length(p1) > 0 & length(p2) > 0) {
        del <- data.frame(p1, p2)
        while (any(del$p1 == del$p2)) {
            del$p2 <- ifelse(del$p1 == del$p2,
                             EpiModel::ssample(which(dat$attr$active == 1 & dat$attr$status != "f"), 1), del$p2)
        }
    }else{
        if(letgo) return(NULL)
    }
    
    ## Discordant edgelist (del)
    del$p1.stat <- dat$attr$status[del$p1]
    del$p2.stat <- dat$attr$status[del$p2]
    serodis <- (del$p1.stat == lab[1] & del$p2.stat == lab[2]) |
        (del$p1.stat == lab[2] & del$p2.stat == lab[1])
    
    return(del[serodis, ])
    
}


#' Infection function
#'
#' Function to implement infection processes for SEIQHRF model
#'
#' @param dat Merged input parameters.
#' @param at Step number
#' @param seed random seed
#'
#' @return Updated dat
#'
#' @importFrom EpiModel ssample
#' @importFrom stats rbinom
#' @importFrom dplyr between
#' @export
infection.FUN <- function(dat, at, seed = NULL){
    
    if(!is.null(seed)) set.seed(seed)
    
    type <- dat$control$type
    nsteps <- dat$control$nsteps
    
    ## Extract param
    param_list <- c("act.rate.i", "inf.prob.i")
    if(type %in% c("SEIQHR", "SEIQHRF")) param_list <- c(param_list, 
                                                         "quar.rate", "disch.rate", "act.rate.e",
                                                         "inf.prob.e", "act.rate.q", "inf.prob.q")
    check_idx <- which(names(dat$param) %in% param_list)
    for(i in check_idx){
        assign(names(dat$param)[i], dat$param[[i]])
    }
    
    # Transmission from infected (a)
    ## Expected acts
    acts <- round(ifelse(length(act.rate.i) > 1, act.rate.i[at - 1], act.rate.i) *
                      dat$epi$num[at - 1] / 2 )
    
    ## Edgelist
    del <- get_del(dat, c("s", "i"), acts)
    
    ## Transmission on edgelist
    if (nrow(del) > 0) {
        
        del$tprob <- ifelse(length(inf.prob.i) > 1, inf.prob.i[at], inf.prob.i)
        
        if (!is.null(dat$param$inter.eff.i) &&
            dplyr::between(at, dat$param$inter.start.i, dat$param$inter.stop.i)) {
            del$tprob <- del$tprob * (1 - dat$param$inter.eff.i)
        }
        
        del$trans <- stats::rbinom(nrow(del), 1, del$tprob)
        del <- del[del$trans == TRUE, ]
        
        if (nrow(del) > 0) {
            newIds <- unique(ifelse(del$p1.stat == "s", del$p1, del$p2)) ## estract individual ID for "s"
            nExp.i <- length(newIds)
            dat$attr$status[newIds] <- "e"
            dat$attr$expTime[newIds] <- at
        } else {
            nExp.i <- 0
        }
        
    } else {
        nExp.i  <- 0
    }
    
    if (type == "SEIQHRF"){
        
        # Transmission from exposed
        ## Expected acts
        acts <- round(ifelse(length(act.rate.e) > 1, act.rate.e[at - 1], act.rate.e) *
                          dat$epi$num[at - 1] / 2 )
        
        ## Edgelist: s, e
        del <- get_del(dat, c("s", "e"), acts, letgo = TRUE)
        if(is.null(del)){
            nExp.e <- 0
        }else{
            ## Transmission on edgelist
            if (nrow(del) > 0) {
                
                del$tprob <- ifelse(length(inf.prob.e) > 1, inf.prob.e[at], inf.prob.e)
                
                if (!is.null(dat$param$inter.eff.e) &&
                    dplyr::between(at, dat$param$inter.start.e, dat$param$inter.stop.e)) {
                    del$tprob <- del$tprob * (1 - dat$param$inter.eff.e)
                }
                
                del$trans <- stats::rbinom(nrow(del), 1, del$tprob)
                del <- del[del$trans == TRUE, ]
                
                if (nrow(del) > 0) {
                    newIds <- unique(ifelse(del$p1.stat == "s", del$p1, del$p2))
                    nExp.e <- length(newIds)
                    dat$attr$status[newIds] <- "e"
                    dat$attr$expTime[newIds] <- at
                } else {
                    nExp.e <- 0
                }
                
            } else {
                nExp.e <- 0
            }
        }
        
        # Transmission from quarantined
        ## Expected acts
        acts <- round(ifelse(length(act.rate.q) > 1, act.rate.q[at - 1], act.rate.q) *
                          dat$epi$num[at - 1] / 2 )
        
        ## Edgelist:s, q
        del <- get_del(dat, c("s", "q"), acts, letgo = TRUE, seed)
        if(is.null(del)){
            nExp.q <- 0
        }else{
            if (nrow(del) > 0) {
                
                del$tprob <- ifelse(length(inf.prob.q) > 1, inf.prob.q[at], inf.prob.q)
                
                if (!is.null(dat$param$inter.eff.q) &&
                    dplyr::between(at, dat$param$inter.start.q, dat$param$inter.stop.q)) {
                    del$tprob <- del$tprob * (1 - dat$param$inter.eff.q)
                }
                
                del$trans <- stats::rbinom(nrow(del), 1, del$tprob)
                del <- del[del$trans == TRUE, ]
                if (nrow(del) > 0) {
                    newIds <- unique(ifelse(del$p1.stat == "s", del$p1, del$p2))
                    nExp.q <- length(newIds)
                    dat$attr$status[newIds] <- "e"
                    dat$attr$expTime[newIds] <- at
                } else {
                    nExp.q <- 0
                }
            } else {
                nExp.q <- 0
            }
        }
        
    }
    
    ## Output
    tmp <- ifelse(type %in% c("SEIQHR", "SEIQHRF"), nExp.i + nExp.q, nExp.i)
    if (at == 2) {
        dat$epi$se.flow <- c(0, tmp)
    } else {
        dat$epi$se.flow[at] <- tmp
    }
    
    dat
}