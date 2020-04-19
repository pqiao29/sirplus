#' Initialize ICM
#'
#' Function to initialize the icm
#'
#' @param param ICM parameters.
#' @param init Initial value parameters.
#' @param control Control parameters
#' @param seed random seed for checking consistency with other versions.
#'
#' @return Updated dat
#' @export
initialize.FUN <- function(param, init, control, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  
  ## Master List for Data ##
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control
  
  # Set attributes
  dat$attr <- list()
  numeric.init <- init[which(sapply(init, class) == "numeric")]
  n <- do.call("sum", numeric.init)
  dat$attr$active <- rep(1, n)
  dat$attr$group <- rep(1, n)
  
  # Initialize status and infection time
  dat <- init_status.icm(dat)
  
  # Summary out list
  dat <- get_prev.icm(dat, at = 1)
  
  return(dat)
}


#' Initialize status of ICM
#'
#' Function to get the status of the initialized icm
#'
#' @param dat Object containing all data
#'
#' @return Updated dat
#' @importFrom EpiModel ssample
#' @export
init_status.icm <- function(dat) {
  
  # Variables ---------------------------------------------------------------
  type <- dat$control$type
  group <- dat$attr$group
  nGroups <- dat$param$groups
  
  nG <- sum(group == 1)
  
  e.num <- dat$init$e.num
  i.num <- dat$init$i.num
  q.num <- dat$init$q.num
  h.num <- dat$init$h.num
  r.num <- dat$init$r.num
  f.num <- dat$init$f.num
  
  # Status ------------------------------------------------------------------
  status <- rep("s", nG)
  status[sample(which(group == 1), size = i.num)] <- "i"
  if (type %in% c("SIR", "SEIR", "SEIQHR", "SEIQHRF")) {
    status[sample(which(group == 1 & status == "s"), size = r.num)] <- "r"
  }
  if (type %in% c("SEIR", "SEIQHR", "SEIQHRF")) {
    status[sample(which(group == 1 & status == "s"), size = e.num)] <- "e"
  }
  if (type %in% c("SEIQHR", "SEIQHRF")) {
    status[sample(which(group == 1 & status == "s"), size = q.num)] <- "q"
    status[sample(which(group == 1 & status == "s"), size = h.num)] <- "h"
  }
  if (type %in% c("SEIQHRF")) {
    status[sample(which(group == 1 & status == "s"), size = f.num)] <- "f"
  }
  
  dat$attr$status <- status
  n <- length(status)
  
  # leave exposure time uninitialised for now, and
  # just set to NA at start.
  
  # Exposure Time ----------------------------------------------------------
  dat$attr$expTime <- rep(NA, n)
  
  # Infection Time ----------------------------------------------------------
  infTime <- rep(NA, n)
  infTime[status == "i"] <- 1
  dat$attr$infTime <- infTime
  
  # Recovery Time ----------------------------------------------------------
  dat$attr$recovTime <- rep(NA, n)
  
  # Need for Hospitalisation Time ----------------------------------------------------------
  dat$attr$hospTime <- rep(NA, n)
  
  # Quarantine Time ----------------------------------------------------------
  dat$attr$quarTime <- rep(NA, n)
  
  # Hospital-need cessation  Time ----------------------------------------------------------
  dat$attr$dischTime <- rep(NA, n)
  
  # Case-fatality  Time ----------------------------------------------------------
  dat$attr$fatTime <- rep(NA, n)
  
  return(dat)
}
