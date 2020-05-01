#' Arrivals function
#'
#' Function to handle background demographics for the SEIQHRF model.
#' Specifically arrivals (births and immigration). Uses the original EpiModel
#' code currently. A replacement that implements modelling the arrival of
#' infected individuals is under development -- but for now, all arrivals
#' go into the S compartment..
#'
#' @param dat Input data needed.
#' @param at time point
#' @param seed random seed for checking consistency
#'
#' @return Updated dat
#' @export
arrivals.FUN <- function(dat, at, seed = NULL) {
  
  if(!is.null(seed)) set.seed(seed)
  # Conditions --------------------------------------------------------------
  if (!dat$param$vital) return(dat)
  
  
  
  flare.inf.point <- dat$param$flare.inf.point
  flare.inf.num <- dat$param$flare.inf.num
  
  if(!(at %in% flare.inf.point)){
    # Variables ---------------------------------------------------------------
    a.rate <- dat$param$a.rate
    a.prop.e <- dat$param$a.prop.e
    a.prop.i <- dat$param$a.prop.i
    a.prop.q <- dat$param$a.prop.q
    a.rand <- dat$control$a.rand
    nOld <- dat$epi$num[at - 1]
    # Process: partition arrivals into compartments -----------------------------------------------------------------
    nArrivals <- ifelse(a.rand, sum(stats::rbinom(nOld, 1, a.rate)), round(nOld * a.rate))
    nArrivals.e <- round(nArrivals * ifelse(length(a.prop.e) > 1, a.prop.e[at], a.prop.e))
    nArrivals.i <- round(nArrivals * ifelse(length(a.prop.i) > 1, a.prop.i[at], a.prop.i))
    nArrivals.q <- round(nArrivals * ifelse(length(a.prop.q) > 1, a.prop.q[at], a.prop.q))
    totArrivals <- nArrivals.e + nArrivals.i + nArrivals.q
    dat$attr$status <- c(dat$attr$status,
                         rep("e", nArrivals.e),
                         rep("i", nArrivals.i),
                         rep("q", nArrivals.q))
    
  }else{
    nArrivals.i <- totArrivals <- flare.inf.num[which(flare.inf.point == at)]
    dat$attr$status <- c(dat$attr$status, rep("i", nArrivals.i) )
    nArrivals.e <- nArrivals.q <- 0
  }
  
  dat$attr$active <- c(dat$attr$active, rep(1, totArrivals))
  dat$attr$group <- c(dat$attr$group, rep(1, totArrivals))
  dat$attr$expTime <- c(dat$attr$expTime, rep(NA, totArrivals))
  dat$attr$infTime <- c(dat$attr$infTime, rep(NA, totArrivals))
  dat$attr$quarTime <- c(dat$attr$quarTime, rep(NA, totArrivals))
  dat$attr$hospTime <- c(dat$attr$ihospTime, rep(NA, totArrivals))
  dat$attr$recovTime <- c(dat$attr$recovTime, rep(NA, totArrivals))
  dat$attr$fatTime <- c(dat$attr$fatTime, rep(NA, totArrivals))
  
  
  # Output ------------------------------------------------------------------
  if (at == 2) {
    dat$epi$a.flow <- c(0, totArrivals)
    dat$epi$a.e.flow <- c(0, nArrivals.e)
    dat$epi$a.i.flow <- c(0, nArrivals.i)
    dat$epi$a.q.flow <- c(0, nArrivals.q)
  } else {
    dat$epi$a.flow[at] <- totArrivals
    dat$epi$a.e.flow[at] <- c(0, nArrivals.e)
    dat$epi$a.i.flow[at] <- c(0, nArrivals.i)
    dat$epi$a.q.flow[at] <- c(0, nArrivals.q)
  }
  
  return(dat)
}