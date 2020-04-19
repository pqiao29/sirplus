#' Progress icm
#'
#' Function to get progress of icms
#'
#' @param dat Object containing all data
#' @param at time point
#' @param seed random seed for checking consistency
#'
#' @return progress
#' @importFrom stats pweibull
#' @importFrom stats rbinom
#' @importFrom EpiModel ssample
#' @importFrom stats rgeom
#' @importFrom stats sd
#' @export
recovery.FUN <- function(dat, at, seed = NULL) {
  
  #cat("at: ", at, "\n")
  
  if(!is.null(seed)) set.seed(seed)
  # Conditions --------------------------------------------------------------------------
  if (!(dat$control$type %in% c("SIR", "SIS", "SEIR", "SEIQHR", "SEIQHRF"))) {
    return(dat)
  }
  
  # Variables ---------------------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status
  groups <- dat$param$groups
  group <- dat$attr$group
  type <- dat$control$type
  expTime <- dat$attr$expTime
  
  # ------------------------- progress from exposed to infectious --------------------------
  res <- update_status(rate = dat$param$prog.rate,
                       rand = dat$control$prog.rand, 
                       active, dat$attr$status, 
                       label = "e", 
                       state =  "i", 
                       at, expTime, 
                       prog.dist.scale = dat$param$prog.dist.scale,
                       prog.dist.shape = dat$param$prog.dist.shape)
  
  nProg <- res[[1]]
  if(!is.null(res[[2]])) dat$attr$infTime[res[[2]]] <- at
  status <- res[[3]]
  dat$attr$status <- status
  
  #cat("finished prog \n")
  
  if (type %in% c("SEIQHR", "SEIQHRF")) {
    
    # ----- quarantine ------- 
    quar.rate <- dat$param$quar.rate
    rate <- ifelse(length(quar.rate) > 1, quar.rate[at], quar.rate)
    
    res <- update_status(rate,
                         rand = dat$control$quar.rand,
                         active = dat$attr$active, 
                         status = dat$attr$status,
                         label = "i", 
                         state = "q",
                         at, expTime, 
                         prog.dist.scale = dat$param$quar.dist.scale,
                         prog.dist.shape = dat$param$quar.dist.shape)
    nQuar <- res[[1]]
    if(!is.null(res[[2]])) dat$attr$quarTime[res[[2]]] <- at
    status <- res[[3]]
    dat$attr$status <- status
    
    # ----- need to be hospitalised -------
    rate <- rep(dat$param$hosp.rate, 
                sum(dat$attr$active == 1 & (dat$attr$status == "i" | dat$attr$status == "q")))
    
    res <- update_status(rate,
                         rand = dat$control$hosp.rand,
                         active = dat$attr$active,
                         status = dat$attr$status, 
                         label = c("i", "q"), 
                         state = "h", 
                         at, expTime, 
                         prog.dist.scale = dat$param$hosp.dist.scale,
                         prog.dist.shape = dat$param$hosp.dist.shape)
    
    nHosp <- res[[1]]
    dat$attr$hospTime[res[[2]]] <- at
    status <- res[[3]]
    dat$attr$status <- status
    
    # ------------------------- discharge from need to be hospitalised ---------------------------
    recovState <- ifelse(type %in% c("SIR", "SEIR", "SEIQHR", "SEIQHRF"), "r", "s")
    disch.rate <- dat$param$disch.rate
    dcrate <- ifelse(length(disch.rate) > 1, disch.rate[at], disch.rate)
    
    res <- update_status(rate = dcrate,
                         rand = dat$control$disch.rand,
                         active = dat$attr$active, 
                         status = dat$attr$status,
                         label = "h", 
                         state = recovState, 
                         at, expTime, 
                         prog.dist.scale = dat$param$disch.dist.scale,
                         prog.dist.shape = dat$param$disch.dist.shape)
    
    nDisch <- res[[1]]
    dat$attr$dischTime[res[[2]]] <- at
    status <- res[[3]]
    dat$attr$status <- status
    #cat("finished discharge \n")
  }
  
  # --------------------------------------- recover -------------------------------------------------
  res <- update_status(rate = dat$param$rec.rate,
                       rand = dat$control$rec.rand, 
                       active, 
                       status = dat$attr$status, 
                       label = c("i", "q"), 
                       state = recovState, 
                       at, 
                       expTime = dat$attr$expTime, 
                       prog.dist.scale = dat$param$rec.dist.scale,
                       prog.dist.shape = dat$param$rec.dist.shape)
  
  nRecov <- res[[1]]
  dat$attr$recovTime[res[[2]]] <- at
  status <- res[[3]]
  dat$attr$status <- status
  
  
  res <- update_status(rate = dat$param$arec.rate,
                       rand = dat$control$arec.rand, 
                       active, 
                       status = dat$attr$status, 
                       label = c("e"), 
                       state = recovState, 
                       at, 
                       expTime = dat$attr$expTime, 
                       prog.dist.scale = dat$param$arec.dist.scale,
                       prog.dist.shape = dat$param$arec.dist.shape)
  
  nRecov <- nRecov + res[[1]]
  dat$attr$recovTime[res[[2]]] <- at
  status <- res[[3]]
  dat$attr$status <- status
  
  if (type %in% c("SEIQHRF")) {
    # -------------------------------------- case fatality ----------------------------------------
    ## obtain rates
    idsElig <- which(active == 1 & status =="h")
    if(length(idsElig) > 0){
      hosp.cap <- dat$param$hosp.cap
      rates <- dat$param$fat.rate.base
      gElig <- group[idsElig]
      timeInHospElig <- at - dat$attr$hospTime[idsElig]
      h.num.yesterday <- 0
      if (!is.null(dat$epi$h.num[at - 1])) {
        h.num.yesterday <- dat$epi$h.num[at - 1]
        if (h.num.yesterday > hosp.cap) {
          blended.rate <- ((hosp.cap * rates) +
                             ((h.num.yesterday - hosp.cap) * dat$param$fat.rate.overcap)) /
            h.num.yesterday
          rates <- blended.rate
        }
      }
      ratesElig <- rates[gElig]
      ratesElig <- ratesElig + timeInHospElig*dat$param$fat.tcoeff*ratesElig
    }
    
    res <- update_status(rate = ratesElig,
                         rand = dat$control$fat.rand,
                         active = dat$attr$active, 
                         status = dat$attr$status,
                         label = "h", 
                         state = "f", at)
    nFat <- res[[1]]
    if(!is.null(res[[2]])) dat$attr$fatTime[res[[2]]] <- at
    status <- res[[3]]
    dat$attr$status <- status
    
    #cat("finished fatality \n")
  }
  
  # Output ------------------------------------------------------------------
  outName_a <- ifelse(type %in% c("SIR", "SEIR"), "ir.flow", "is.flow")
  if (type %in% c("SEIR", "SEIQHR", "SEIQHRF")) {
    outName_b <- "ei.flow"
  }
  if (type %in% c("SEIQHR", "SEIQHRF")) {
    outName_c <- "iq.flow"
    outName_d <- "iq2h.flow"
  }
  if (type %in% c("SEIQHRF")) {
    outName_e <- "hf.flow"
  }
  ## Summary statistics
  if (at == 2) {
    dat$epi[[outName_a[1]]] <- c(0, nRecov)
    if (type %in% c("SEIR", "SEIQHR")) {
      dat$epi[[outName_b[1]]] <- c(0, nProg)
    }
    if (type %in% c("SEIQHR", "SEIQHRF")) {
      dat$epi[[outName_c[1]]] <- c(0, nQuar)
      dat$epi[[outName_d[1]]] <- c(0, nHosp)
    }
    if (type %in% c("SEIQHRF")) {
      dat$epi[[outName_e[1]]] <- c(0, nFat)
    }
  } else {
    dat$epi[[outName_a[1]]][at] <- nRecov
    if (type %in% c("SEIR", "SEIQHR")) {
      dat$epi[[outName_b[1]]][at] <- nProg
    }
    if (type %in% c("SEIQHR", "SEIQHRF")) {
      dat$epi[[outName_c[1]]][at] <- nQuar
      dat$epi[[outName_d[1]]][at] <- nHosp
    }
    if (type %in% c("SEIQHRF")) {
      dat$epi[[outName_e[1]]][at] <- nFat
    }
  }
  
  return(dat)
}
