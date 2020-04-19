#' Control Settings for Stochastic Individual Contact Models
#'
#' Sets the controls for stochastic individual contact models simulated 
#' with \code{\link{simulate_seiqhrf}}. Similar to EpiModel::control.icm, but allows for
#' model types with additional compartments (e.g. 'SEIQHRF').
#'
#' @param type Disease type to be modeled, with the choice of \code{"SI"} for
#'        Susceptible-Infected diseases, \code{"SIR"} for
#'        Susceptible-Infected-Recovered diseases, and \code{"SIS"} for
#'        Susceptible-Infected-Susceptible diseases.
#' @param nsteps Number of time steps to solve the model over. This must be a
#'        positive integer.
#' @param nsims Number of simulations to run.
#' @param prog.rand Method for progression from E compartment to I. If TRUE, 
#'         random binomial draws at prog.rate, if FALSE, random draws from a 
#'         Weibull distribution (yes, I know it should be a discrete Weibull 
#'         distribution but it makes little difference and speed of computation 
#'         matters), with parameters prog.dist.scale and prog.dist.shape
#' @param rec.rand If \code{TRUE}, use a stochastic recovery model, with the
#'        number of recovered at each time step a function of random draws from
#'        a binomial distribution with the probability equal to \code{rec.rate}.
#'        If \code{FALSE}, then a deterministic rounded count of the expectation
#'        implied by that rate.
#' @param quar.rand Method for self-isolation transition from I to Q. If TRUE, 
#'          random binomial draws at quar.rate, if FALSE, random sample with a 
#'          sample fraction also given by `quar.rate.
#' @param hosp.rand Method for transition from I or Q to H -- that is, from 
#'          infectious or from self-isolated to requiring hospitalisation. If 
#'          TRUE, random binomial draws at hosp.rate, if FALSE, random sample 
#'          with a sample fraction also given by `hosp.rate.
#' @param disch.rand Method for transition from H to R -- that is, from 
#'          requiring hospitalisation to recovered. If TRUE, random binomial 
#'          draws at disch.rate, if FALSE, random sample with a sample fraction 
#'          also given by disch.rate. Note that the only way out of the H 
#'          compartment is recovery or death.
#' @param fat.rand Method for case fatality transition from H to F. If TRUE, 
#'          random binomial draws at fat.rate.base, if FALSE, random sample with
#'          a sample fraction also given by fat.rate.base. However, if the 
#'          current number of patients in the H (needs hospitalisation) 
#'          compartment is above a hospital capacity level specified by 
#'          hosp.cap, then the fatality rate is the mean of the base fatality 
#'          rate weighted by the hospital capacity, plus a higher rate, 
#'          specified by fat.rate.overcap, weighted by the balance of those 
#'          requiring hospitalisation (but who can't be accommodated). By 
#'          setting fat.rate.overcap higher, the effect of exceeding the 
#'          capacity of the health care system can be simulated. There is also
#'          a coefficient fat.tcoeff for the fatality rates that increases them 
#'          as a linear function of the number of days the individual has been 
#'          in the H compartment. Use of the co-efficient better approximates 
#'          the trapezoid survival time distribution typical of ICU patients.
#' @param arec.rand  Method for recovery transition from E to R. If TRUE, 
#'         random binomial draws at arec.rate, if FALSE, random draws from a 
#'         random draws from a Weibull distribution, with parameters 
#'         arec.dist.scale and arec.dist.shape.
#' @param a.rand If \code{TRUE}, use a stochastic arrival model, with the
#'        number of arrivals at each time step a function of random draws from a
#'        binomial distribution with the probability equal to the governing arrival
#'        rates. If \code{FALSE}, then a deterministic rounded count of the
#'        expectation implied by those rates.
#' @param d.rand If \code{TRUE}, use a stochastic departure model, with the number of
#'        departures at each time step a function of random draws from a binomial
#'        distribution with the probability equal to the governing departure rates.
#'        If \code{FALSE}, then a deterministic rounded count of the expectation
#'        implied by those rates.
#' @param initialize.FUN Module to initialize the model at the outset, with the
#'        default function of \code{\link{initialize.icm}}.
#' @param infection.FUN Module to simulate disease infection, with the default
#'        function of \code{\link{infection.icm}}.
#' @param recovery.FUN Module to simulate disease recovery, with the default
#'        function of \code{\link{recovery.icm}}.
#' @param departures.FUN Module to simulate departures or exits, with the default
#'        function of \code{\link{departures.icm}}.
#' @param arrivals.FUN Module to simulate arrivals or entries, with the default
#'        function of \code{\link{arrivals.icm}}.
#' @param get_prev.FUN Module to calculate disease prevalence at each time step,
#'        with the default function of \code{\link{get_prev.icm}}.
#' @param verbose If \code{TRUE}, print model progress to the console.
#' @param verbose.int Time step interval for printing progress to console, where
#'        0 (the default) prints completion status of entire simulation and
#'        positive integer \code{x} prints progress after each \code{x} time
#'        steps.
#' @param skip.check If \code{TRUE}, skips the default error checking for the
#'        structure and consistency of the parameter values, initial conditions,
#'        and control settings before running base epidemic models. Setting
#'        this to \code{FALSE} is recommended when running models with new modules
#'        specified.
#' @param ncores Number of physical CPU cores used for parallel computation.
#' @param ... Additional control settings passed to model.
#'
#' @details
#' \code{control} sets the required control settings for any stochastic
#' individual contact model solved with the \code{simulate_seiqhrf} function. Controls
#' are required for both base model types and when passing original process
#' modules. For an overview of control settings for base ICM class models,
#' consult the \href{http://statnet.github.io/tut/BasicICMs.html}{Basic ICMs}
#' tutorial. For all base models, the \code{type} argument is a necessary
#' parameter and it has no default.
#'
#' @section New Modules:
#' Base ICM models use a set of module functions that specify
#' how the individual agents in the population are subjected to infection, recovery,
#' demographics, and other processes. Core modules are those listed in the
#' \code{.FUN} arguments. For each module, there is a default function used in
#' the simulation. The default infection module, for example, is contained in
#' the \code{\link{infection.FUN}} function.
#'
#' For original models, one may substitute replacement module functions for any of
#' the default functions. New modules may be added to the workflow at each time
#' step by passing a module function via the \code{...} argument.
#'
#' @keywords parameterization
#'
#' @return A list of control parameters and core functions
#' @export
control_seiqhrf <- function(type = "SEIQHRF", nsteps = 366, nsims = 8, 
                            prog.rand = FALSE, quar.rand = TRUE, hosp.rand = TRUE,
                            disch.rand = TRUE, rec.rand = FALSE, arec.rand = TRUE,
                            fat.rand = TRUE, a.rand = TRUE, d.rand = TRUE, 
                            initialize.FUN = 'initialize.FUN', infection.FUN = 'infection.FUN', 
                            recovery.FUN = 'recovery.FUN', departures.FUN = 'departures.FUN',
                            arrivals.FUN = 'arrivals.FUN', get_prev.FUN = 'get_prev.FUN',
                            verbose = FALSE, verbose.int = 0, skip.check = FALSE, 
                            ncores = 4, ...) {

  # Get arguments
  p <- list()
  passed <- names(as.list(match.call())[-1])
  p$usr.specified <- passed
  
  formal.args <- formals(sys.function())
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    if (as.logical(mget(arg) != "")) {
      p[arg] <- list(get(arg))
    }
  }
  dot.args <- list(...)
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }

  if ("births.FUN" %in% names(dot.args)) {
    p$arrivals.FUN <- dot.args$births.FUN
    p$births.FUN <- dot.args$births.FUN <- NULL
    message("EpiModel 1.7.0 onward renamed the birth function births.FUN to 
            arrivals.FUN. See documentation for details.")
  }
  if ("deaths.FUN" %in% names(dot.args)) {
    p$departures.FUN <- dot.args$deaths.FUN
    p$deaths.FUN <- dot.args$deaths.FUN <- NULL
    message("EpiModel 1.7.0 onward renamed the death function deaths.FUN to 
            departures.FUN. See documentation for details.")
  }


  ## Module classification
  p$bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  #p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)


  ## Defaults and checks
  if (is.null(p$type) | !(p$type %in% c("SI", "SIS", "SIR", "SEIR", "SEIQHR", 
                                        "SEIQHRF"))) {
    stop("Specify type as \"SI\", \"SIS\", \"SIR\", \"SEIR\", \"SEIQHR\", or 
         \"SEIQHRF\" ", call. = FALSE)
  }
  if (is.null(p$nsteps)) stop("Specify nsteps", call. = FALSE)
  if (is.null(p$nsims)) stop("Specify nsims", call. = FALSE)

  ## Output
  class(p) <- c("control.seiqhrf", "list")
  return(p)
}
