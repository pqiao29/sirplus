#' Merge icm
#'
#' Function to merge icms
#'
#' @param x First to merge
#' @param y Second to merge
#' @param ... other parameters
#'
#' @return Combined
#' 
#' @importFrom doParallel registerDoParallel
#' 
merge.seiqhrf.icm <- function(x, y, ...) {

  ## Check structure
  if (length(x) != length(y) || !identical(names(x), names(y))) {
    stop("x and y have different structure")
  }
  if (x$control$nsims > 1 & y$control$nsims > 1 &
      !identical(sapply(x, class), sapply(y, class))) {
    stop("x and y have different structure")
  }

  ## Check params
  check1 <- identical(x$param, y$param)
  check2 <- identical(x$control[-which(names(x$control) %in% c("nsims", "currsim"))],
                      y$control[-which(names(y$control) %in% c("nsims", "currsim"))])

  if (check1 == FALSE) {
    stop("x and y have different parameters")
  }
  if (check2 == FALSE) {
    stop("x and y have different controls")
  }


  z <- x
  new.range <- (x$control$nsims + 1):(x$control$nsims + y$control$nsims)

  # Merge data
  for (i in 1:length(x$epi)) {
    if (x$control$nsims == 1) {
      x$epi[[i]] <- data.frame(x$epi[[i]])
    }
    if (y$control$nsims == 1) {
      y$epi[[i]] <- data.frame(y$epi[[i]])
    }
    z$epi[[i]] <- cbind(x$epi[[i]], y$epi[[i]])
    names(z$epi[[i]])[new.range] <- paste0("sim", new.range)
  }

  z$control$nsims <- max(new.range)

  return(z)
}

#' SEIQHRF function
#'
#' Function to implement infection processes for SEIQHRF model
#'
#' @param init Initial state values, returned by \code{\link{init_seiqhrf}}.
#' @param control Control parameters of the model, returned by \code{\link{control_seiqhrf}}.
#' @param param Model parameters, returned by \code{\link{param_seiqhrf}}. 
#' @param seed random seed for checking consistency
#'
#' @return Updated dat
#'
#' @importFrom foreach "%dopar%" foreach
#' @importFrom EpiModel verbose.icm
#' @importFrom future availableCores
#' @export
seiqhrf <- function(init, control, param, seed = NULL) {


  crosscheck_seiqhrf(param, init, control)
  nsims <- control$nsims
  ncores <- ifelse(control$nsims == 1, 1, min(future::availableCores(), control$ncores))
  control$ncores <- ncores

  if (ncores == 1) {

      # Simulation loop start
      for (s in 1:control$nsims) {

        ## Initialization module
        if (!is.null(control[["initialize.FUN"]])) {
          dat <- do.call(control[["initialize.FUN"]], list(param, init, control, seed))
        }


        # Timestep loop
        for (at in 2:control$nsteps) {

          ## Infection
          if (!is.null(control[["infection.FUN"]])) {
            dat <- do.call(control[["infection.FUN"]], list(dat, at, seed))
          }


          ## Recovery
          if (!is.null(control[["recovery.FUN"]])) {
            dat <- do.call(control[["recovery.FUN"]], list(dat, at, seed))
          }


          ## Departure Module
          if (!is.null(control[["departures.FUN"]])) {
            dat <- do.call(control[["departures.FUN"]], list(dat, at, seed))
          }


          ## Arrival Module
          if (!is.null(control[["arrivals.FUN"]])) {
            dat <- do.call(control[["arrivals.FUN"]], list(dat, at, seed))
          }


          ## Outputs
          if (!is.null(control[["get_prev.FUN"]])) {
            dat <- do.call(control[["get_prev.FUN"]], list(dat, at))
          }


          ## Track progress
          EpiModel::verbose.icm(dat, type = "progress", s, at)
        }

        # Set output
        if (s == 1) {
          out <- saveout_seiqhrf(dat, s)
        } else {
          out <- saveout_seiqhrf(dat, s, out)
        }

      } # Simulation loop end

  } # end of single core execution

  if (ncores > 1) {
    doParallel::registerDoParallel(ncores)

    sout <- foreach(s = 1:nsims) %dopar% {

        control$nsims <- 1
        control$currsim <- s

        ## Initialization module
        if (!is.null(control[["initialize.FUN"]])) {
          dat <- do.call(control[["initialize.FUN"]], list(param, init, control, seed))
        }

        # Timestep loop
        for (at in 2:control$nsteps) {

          ## Infection
          if (!is.null(control[["infection.FUN"]])) {
            dat <- do.call(control[["infection.FUN"]], list(dat, at, seed))
          }


          ## Recovery
          if (!is.null(control[["recovery.FUN"]])) {
            dat <- do.call(control[["recovery.FUN"]], list(dat, at, seed))
          }


          ## Departure Module
          if (!is.null(control[["departures.FUN"]])) {
            dat <- do.call(control[["departures.FUN"]], list(dat, at, seed))
          }


          ## Arrival Module
          if (!is.null(control[["arrivals.FUN"]])) {
            dat <- do.call(control[["arrivals.FUN"]], list(dat, at, seed))
          }


          ## Outputs
          if (!is.null(control[["get_prev.FUN"]])) {
            dat <- do.call(control[["get_prev.FUN"]], list(dat, at))
          }


          ## Track progress
          EpiModel::verbose.icm(dat, type = "progress", s, at)
        }

        # Set output
        out <- saveout_seiqhrf(dat, s=1)
        class(out) <- "icm"
        return(out)

  }

    # aggregate results collected from each thread
    collected_times <- list()

    # collect the times from sout then delete them
    for (i in 1:length(sout)) {
        collected_times[[paste0("sim", i)]] <- sout[[i]]$times$sim1
        sout[[i]]$times <- NULL
    }

    # merge $epi structures
    merged.out <- sout[[1]]
    for (i in 2:length(sout)) {
      merged.out <- merge.seiqhrf.icm(merged.out, sout[[i]], param.error = FALSE)
    }
    out <- merged.out

    # add the collected timing data
    out$times <- collected_times

  } # end of parallel execution
  
  class(out) <- "seiqhrf"

  return(out)
}