#' Departures function
#'
#' Function to handle background demographics for the SEIQHRF model.
#' Specifically departures (deaths not due to the virus, and emigration).
#'
#' @param dat Merged input parameters.
#' @param at Step number
#' @param seed random seed for checking consistency
#'
#' @return Updated dat
#' @export
#' @examples
#' 
#' control <- control_seiqhrf()
#' param <- param_seiqhrf()
#' init <- init_seiqhrf()
#' 
#' dat <- initialize.FUN(param, init, control) 
#' dat <- infection.FUN(dat, at = 2)  
#' dat <- recovery.FUN(dat, at = 2) 
#' dat <- departures.FUN(dat, at = 2)
#' 

departures.FUN <- function(dat, at, seed = NULL) {
    
    if(!is.null(seed)) set.seed(seed)
    
    # Conditions --------------------------------------------------------------
    if (!dat$param$vital) return(dat)
    
    # Variables -----------------------------------------------------------------
    rate <- dat$param$ds.rate
    rand <- dat$control$d.rand
    status <- dat$attr$status
    active <- dat$attr$active
    
    # Susceptible departures ------------------------------------------------------
    res <- update_active(rate, rand, active, status, label = "s")
    nDepartures <- res[[1]]
    if(!is.null(res[[2]]))  active <- dat$attr$active[res[[2]]] <- 0
    
    if (at == 2) dat$epi$ds.flow <- c(0, nDepartures) else dat$epi$ds.flow[at] <- nDepartures
    
    # Exposed Departures ---------------------------------------------------------
    res <- update_active(rate, rand, active, status, label = "e")
    nDepartures <- res[[1]]
    if(!is.null(res[[2]]))  active <- dat$attr$active[res[[2]]] <- 0
    
    if (at == 2) dat$epi$de.flow <- c(0, nDepartures) else dat$epi$de.flow[at] <- nDepartures
    
    # Infected Departures ---------------------------------------------------------
    res <- update_active(rate, rand, active, status, label = "i")
    nDepartures <- res[[1]]
    if(!is.null(res[[2]]))  active <- dat$attr$active[res[[2]]] <- 0
    
    if (at == 2) dat$epi$di.flow <- c(0, nDepartures) else dat$epi$di.flow[at] <- nDepartures
    
    # Quarantined Departures ---------------------------------------------------------
    res <- update_active(rate, rand, active, status, label = "q")
    nDepartures <- res[[1]]
    if(!is.null(res[[2]]))  active <- dat$attr$active[res[[2]]] <- 0
    
    if (at == 2) dat$epi$dq.flow <- c(0, nDepartures) else dat$epi$dq.flow[at] <- nDepartures
    
    # Hospitalised Departures ---------------------------------------------------------
    res <- update_active(rate, rand, active, status, label = "h")
    nDepartures <- res[[1]]
    if(!is.null(res[[2]]))  active <- dat$attr$active[res[[2]]] <- 0
    
    if (at == 2) dat$epi$dh.flow <- c(0, nDepartures) else dat$epi$dh.flow[at] <- nDepartures
    
    # Recovered Departures --------------------------------------------------------
    res <- update_active(rate, rand, active, status, label = "r")
    nDepartures <- res[[1]]
    if(!is.null(res[[2]]))  active <- dat$attr$active[res[[2]]] <- 0
    
    if (at == 2) dat$epi$dr.flow <- c(0, nDepartures) else dat$epi$dr.flow[at] <- nDepartures
    
    # return --------------------------------------------------------
    return(dat)
}
