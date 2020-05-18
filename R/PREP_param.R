#' @title Epidemic Parameters for Stochastic Individual Contact Models
#'
#' @description Sets the epidemic parameters for stochastic individual contact
#'              models simulated with \code{seiqhrf}.
#' @param inf.prob.e Probability of passing on infection at each
#'     exposure event for interactions between infectious people in the E
#'     compartment and susceptible in S. Note the default is lower than for
#'     inf.prob.i reflecting the reduced infectivity of infected but
#'     asymptomatic people (the E compartment). Otherwise as for inf.exp.i.
#' @param act.rate.e The number of exposure events (acts) between
#'     infectious individuals in the E compartment and susceptible individuals
#'     in the S compartment, per day. Otherwise as for act.rate.i.
#' @param  inf.prob.i Probability of passing on infection at each
#'     exposure event for interactions between infectious people in the I
#'     compartment and susceptible in S. Reducing inf.prob.i is equivalent to
#'     increasing hygiene measures, such as not putting hands in eyes, nose or
#'     moth, use of hand sanitisers, wearing masks by the infected, and so on.
#' @param act.rate.i The number of exposure events (acts) between
#'     infectious individuals in the I compartment and susceptible individuals
#'     in the S compartment, per day. It's stochastic, so the rate is an
#'     average, some individuals may have more or less. Note that not every
#'     exposure event results in infection - that is governed by the inf.prob.i
#'     parameters (see below). Reducing act.rate.i is equivalent to increasing
#'     social distancing by people in the I compartment.
#' @param inf.prob.q  Probability of passing on infection at each
#'     exposure event for interactions between infectious people in the Q
#'     compartment and susceptible in S. Note the default is lower than for
#'     inf.prob.i reflecting the greater care that self-isolated individuals
#'     will, on average, take regarding hygiene measures, such as wearing masks,
#'     to limit spread to others. Otherwise as for inf.exp.i.
#' @param act.rate.q The number of exposure events (acts) between
#'     infectious individuals in the Q compartment (isolated, self or otherwise)
#'     and susceptible individuals in the S compartment, per day. Note the much
#'     lower rate than for the I and E compartments, reflecting the much
#'     greater degree of social isolation for someone in (self-)isolation. The
#'     exposure event rate is not zero for this group, just much less.
#'     Otherwise as for act.rate.i.
#' @param quar.rate Rate per day at which symptomatic (or tested
#'     positive), infected I compartment people enter self-isolation (Q
#'     compartment). Asymptomatic E compartment people can't enter
#'     self-isolation because they don't yet know they are infected. Default is
#'     a low rate reflecting low community awareness or compliance with
#'     self-isolation requirements or practices, but this can be tweaked when
#'     exploring scenarios.
#' @param quar.dist.scale Scale parameter for Weibull distribution for
#'     recovery, see quar.rand for details.
#' @param quar.dist.shape Shape parameter for Weibull distribution for
#'     recovery, see quar.rand for details. Read up on the Weibull distribution
#'     before changing the default.
#' @param hosp.rate Rate per day at which symptomatic (or tested
#'     positive), infected I compartment people or self-isolated Q compartment
#'     people enter the state of requiring hospital care -- that is, become
#'     serious cases. A default rate of 1% per day with an average illness
#'     duration of about 10 days means a bit less than 10% of cases will
#'     require hospitalisation, which seems about right (but can be tweaked,
#'     of course).
#' @param hosp.dist.scale Scale parameter for Weibull distribution for
#'     recovery, see hosp.rand for details.
#' @param hosp.dist.shape Shape parameter for Weibull distribution for
#'     recovery, see hosp.rand for details. Read up on the Weibull distribution
#'     before changing the default.
#' @param disch.rate Rate per day at which people needing
#'     hospitalisation recover.
#' @param disch.dist.scale Scale parameter for Weibull distribution for
#'     recovery, see disch.rand for details.
#' @param disch.dist.shape Shape parameter for Weibull distribution for
#'     recovery, see disch.rand for details. Read up on the Weibull distribution
#'     before changing the default.
#' @param prog.rate Rate per day at which people who are infected
#'     but asymptomatic (E compartment) progress to becoming symptomatic (or
#'     test-positive), the I compartment. See prog.rand above for more details.
#' @param prog.dist.scale Scale parameter for Weibull distribution
#'     for progression, see prog.rand for details.
#' @param prog.dist.shape Shape parameter for Weibull distribution
#'     for progression, see prog.rand for details. Read up on the Weibull
#'     distribution before changing the default.
#' @param rec.rate Rate per day at which people who are infected and
#'     symptomatic (I compartment) recover, thus entering the R compartment.
#'     See rec.rand above for more details.
#' @param rec.dist.scale Scale parameter for Weibull distribution for
#'     recovery, see rec.rand for details.
#' @param rec.dist.shape Shape parameter for Weibull distribution for
#'     recovery, see rec.rand for details. Read up on the Weibull distribution
#'     before changing the default.
#' @param arec.rate Rate per day at which people who are exposed but asymptotic 
#'     (E compartment) recover, thus entering the R compartment.
#'     See arec.rand above for more details.
#' @param arec.dist.scale Scale parameter for Weibull distribution for
#'     recovery, see arec.rand for details.
#' @param arec.dist.shape Shape parameter for Weibull distribution for
#'     recovery, see arec.rand for details. 
#' @param fat.rate.base Baseline mortality rate per day for people
#'     needing hospitalisation (deaths due to the virus). See fat.rand for more
#'     details.
#' @param hosp.cap Number of available hospital beds for the modelled
#'     population. See fat.rand for more details.
#' @param fat.rate.overcap Mortality rate per day for people needing
#'     hospitalisation but who can't get into hospital due to the hospitals
#'     being full (see hosp.cap and fat.rand). The default rate is twice that
#'     for those who do get into hospital.
#' @param fat.tcoeff Time co-efficient for increasing mortality rate
#'     as time in the H compartment increases for each individual in it. See
#'     fat.rand for details.
#' @param vital Enables demographics, that is, arrivals and
#'     departures, to and from the simulated population.
#' @param flare.inf.point (Either a vector or a scalar) Time points where there is a sudden arrival of I's. 
#' @param flare.inf.num (same dimension as flare.inf.point) The number of I's arriving at the specified time points in flare.inf.point.
#' @param a.rate Background demographic arrival rate. Currently all
#'     arrivals go into the S compartment, the default is approximately the
#'     daily birth rate for Australia. Will be extended to cover immigration in
#'     future versions.
#' @param a.prop.e Arrivals proportion that goes to E (immigration). 
#' @param a.prop.i Arrivals proportion that goes to I (immigration). 
#' @param a.prop.q Arrivals proportion that goes to Q (immigration). 
#' @param ds.rate Background demographic departure (death not due to
#'     virus) rates. Defaults based on Australian crude death rates. Can be
#'     used to model emigration as well as deaths.
#' @param de.rate Background demographic departure (death not due to
#'     virus) rates. Defaults based on Australian crude death rates. Can be
#'     used to model emigration as well as deaths.
#' @param di.rate Background demographic departure (death not due to
#'     virus) rates. Defaults based on Australian crude death rates. Can be used
#'     to model emigration as well as deaths.
#' @param dq.rate Background demographic departure (death not due to
#'     virus) rates. Defaults based on Australian crude death rates. Can be used
#'     to model emigration as well as deaths.
#' @param dh.rate Background demographic departure (death not due to
#'     virus) rates. Defaults based on Australian crude death rates. Can be used
#'     to model emigration as well as deaths.
#' @param dr.rate Background demographic departure (death not due to
#'     virus) rates. Defaults based on Australian crude death rates. Can be used
#'     to model emigration as well as deaths.
#' @param ... Other parameters.
#'
#'
#' @details
#' \code{param_seiqhrf} sets the epidemic parameters for the stochastic individual
#' contact models simulated with the \code{\link{seiqhrf}} function. Models
#' may use the base types, for which these parameters are used, or new process
#' modules which may use these parameters (but not necessarily). A detailed
#' description of ICM parameterization for base models is found in the
#' \href{http://statnet.github.io/tut/BasicICMs.html}{Basic ICMs} tutorial.
#'
#' For base models, the model specification will be chosen as a result of
#' the model parameters entered here and the control settings in
#' \code{\link{control_seiqhrf}}.
#' 
#' @section New Modules:
#' To build original models outside of the base models, new process modules
#' may be constructed to replace the existing modules or to supplement the existing
#' set. These are passed into the control settings in \code{\link{control_seiqhrf}}.
#' New modules may use either the existing model parameters named here, an
#' original set of parameters, or a combination of both. The \code{...} allows
#' the user to pass an arbitrary set of original model parameters into
#' \code{param_seiqhrf}. Whereas there are strict checks with default modules for
#' parameter validity, these checks are the user's responsibility with new modules.
#'
#' @seealso Use \code{\link{init.icm}} to specify the initial conditions and
#'          \code{\link{control_seiqhrf}} to specify the control settings. Run the
#'          parameterized model with \code{\link{seiqhrf}}.
#'
#' @keywords parameterization
#'
#' @export
#'
param_seiqhrf <- function(inf.prob.e = 0.02,
                          act.rate.e = 10,
                          inf.prob.i = 0.05,
                          act.rate.i = 10,
                          inf.prob.q = 0.02,
                          act.rate.q = 2.5,
                          prog.rate = 1/10,
                          quar.rate = 1/30,
                          hosp.rate = 1/100,
                          disch.rate = 1/15,
                          rec.rate = 0.071,  
                          arec.rate = 0.05,
                          prog.dist.scale = 5,
                          prog.dist.shape = 1.5,
                          quar.dist.scale = 1,
                          quar.dist.shape = 1,
                          hosp.dist.scale = 1,
                          hosp.dist.shape = 1,
                          disch.dist.scale = 1,
                          disch.dist.shape = 1,
                          rec.dist.scale = 35,
                          rec.dist.shape = 1.5,
                          arec.dist.scale = 35,
                          arec.dist.shape = 1.5,
                          fat.rate.base = 1/50,
                          hosp.cap = 40,
                          fat.rate.overcap = 1/25,
                          fat.tcoeff = 0.5,
                          vital = TRUE,
                          flare.inf.point = NULL, 
                          flare.inf.num = NULL, 
                          a.rate = (10.5/365)/1000,
                          a.prop.e = 0.01,
                          a.prop.i = 0.001,
                          a.prop.q = 0.01,
                          ds.rate = (7/365)/1000,
                          de.rate = (7/365)/1000,
                          di.rate = (7/365)/1000,
                          dq.rate = (7/365)/1000,
                          dh.rate = (20/365)/1000,
                          dr.rate = (7/365)/1000,  ...) {
    
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
    
    if ("b.rate" %in% names.dot.args) {
        p$a.rate <- dot.args$b.rate
        message("EpiModel 1.7.0 onward renamed the birth rate parameter b.rate to a.rate. ",
                "See documentation for details.")
    }
    
    if ("b.rand" %in% names.dot.args) {
        p$a.rand <- dot.args$b.rand
        message("EpiModel 1.7.0 onward renamed the stochastic birth flag b.rand to a.rand. ",
                "See documentation for details.")
    }
    
    ## Defaults and checks
    if (is.null(p$act.rate)) {
        p$act.rate <- 1
    }
    p$vital <- ifelse(!is.null(p$a.rate) | !is.null(p$ds.rate) |
                          !is.null(p$di.rate) | !is.null(p$dr.rate), TRUE, FALSE)
    
    p$groups <- ifelse(any(grepl(".g2", names(p))) == TRUE, 2, 1)
    
    if (!is.null(p$inter.eff) && is.null(p$inter.start)) {
        p$inter.start <- 1
    }
    
    ## Output
    class(p) <- c("param.seiqhrf", "list")
    return(p)
}


