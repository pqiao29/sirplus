#' SEIQHRF Simulation Wrapper
#'
#' Wrapper to implement SEIQHRF model
#'
#' @param type Type of model: SI, SIR, SIS, SEIR, SEIQHR and
#'     SEIQHRF available, but only SEIQHRF is likely to work in the
#'     current version of the code.
#' @param nsteps Number of time steps to solve the model over. This must be a
#'         positive integer.
#' @param nsims Number of simulations to run.
#' @param ncores Number of physical CPU cores used for parallel computation.
#' @param prog.rand Method for progression from E compartment to I. If TRUE, 
#'         random binomial draws at prog.rate, if FALSE, random draws from a 
#'         Weibull distribution (yes, I know it should be a discrete Weibull 
#'         distribution but it makes little difference and speed of computation 
#'         matters), with parameters prog.dist.scale and prog.dist.shape
#' @param rec.rand  Method for recovery transition from I, Q or H to R. If TRUE, 
#'         random binomial draws at rec.rate, if FALSE, random draws from a 
#'         random draws from a Weibull distribution, with parameters 
#'         rec.dist.scale and rec.dist.shape.
#' @param arec.rand  Method for recovery transition from E to R. If TRUE, 
#'         random binomial draws at arec.rate, if FALSE, random draws from a 
#'         random draws from a Weibull distribution, with parameters 
#'         arec.dist.scale and arec.dist.shape.
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
#' @param infection.FUN No, being infected with SARS-CoV2 is not fun. Rather 
#'          this is the the name of the function to implement infection 
#'          processes. Use the default.
#' @param recovery.FUN Function to implement recovery processes. Use the 
#'          default.
#' @param departures.FUN Handles background demographics, specifically 
#'           departures (deaths not due to the virus, and emigration). Use the 
#'           default.
#' @param arrivals.FUN Handles background demographics, specifically arrivals 
#'           (births and immigration). Uses the original EpiModel code 
#'           currently. A replacement that implements modelling the arrival of 
#'           infected individuals is under development -- but for now, all 
#'           arrivals go into the S compartment.
#' @param get_prev.FUN Utility function that collects prevalence and transition 
#'          time data from each run and stores it away in the simulation result 
#'          object. Use the default.
#' @param s.num Initial number of *S compartment individuals in
#'     the simulated population. An overall population of 10,000 is a good
#'     compromise. A set of models will still take several minutes or more
#'     to run, in parallel.
#' @param e.num Initial number of E compartment individuals in
#'     the simulated population.
#' @param i.num Initial number of I compartment individuals in
#'     the simulated population.
#' @param q.num Initial number of Q compartment individuals in
#'     the simulated population.
#' @param h.num Initial number of H compartment individuals in
#'      the simulated population.
#' @param r.num Initial number of R compartment individuals in
#'     the simulated population.
#' @param f.num Initial number of F compartment individuals in
#'     the simulated population.
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
#' @param out Summary function for the simulation runs. median is
#'     also available, or percentiles, see the EpiModel documentation.
#' @param seed Random see for checking consistency.    
#'
#' @return list with simulation and simulation dataframe
#' @export
simulate_seiqhrf <- function(type = "SEIQHRF",
    nsteps = 366,
    nsims = 8,
    ncores = 4,
    prog.rand = FALSE,
    quar.rand = TRUE,
    hosp.rand = TRUE,
    disch.rand = TRUE,
    rec.rand = FALSE,
    arec.rand = TRUE,
    fat.rand = TRUE,
    infection.FUN = 'infection.FUN',  # infection.seiqhrf.icm,
    recovery.FUN = 'recovery.FUN', # progress.seiqhrf.icm,
    departures.FUN = 'departures.FUN', # departures.seiqhrf.icm,
    arrivals.FUN = 'arrivals.FUN', # arrivals.seiqhrf.icm,
    get_prev.FUN =  'get_prev.FUN', # get_prev.seiqhrf.icm,
    # init.icm params
    s.num = 9997,
    e.num=0,
    i.num = 3,
    q.num=0,
    h.num=0,
    r.num = 0,
    f.num = 0,
    # param.icm params
    inf.prob.e = 0.02,
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
    a.rate = (10.5/365)/1000,
    a.prop.e = 0.01,
    a.prop.i = 0.001,
    a.prop.q = 0.01,
    ds.rate = (7/365)/1000,
    de.rate = (7/365)/1000,
    di.rate = (7/365)/1000,
    dq.rate = (7/365)/1000,
    dh.rate = (20/365)/1000,
    dr.rate = (7/365)/1000,
    out="mean", seed = NULL
) {

    control <- control_seiqhrf(type = type,
                               nsteps = nsteps,
                               nsims = nsims,
                               ncores = ncores,
                               prog.rand = prog.rand,
                               quar.rand = quar.rand, 
                               hosp.rand = hosp.rand, 
                               disch.rand = disch.rand, 
                               rec.rand = rec.rand,
                               arec.rand = arec.rand, 
                               infection.FUN = infection.FUN,
                               recovery.FUN = recovery.FUN,
                               arrivals.FUN = arrivals.FUN,
                               departures.FUN = departures.FUN,
                               get_prev.FUN = get_prev.FUN)

    init <- init_seiqhrf(s.num = s.num,
                     e.num = e.num,
                     i.num = i.num,
                     q.num = q.num,
                     h.num = h.num,
                     r.num = r.num,
                     f.num = f.num)

    param <- param_seiqhrf(inf.prob.e = inf.prob.e,
                        act.rate.e = act.rate.e,
                        inf.prob.i = inf.prob.i,
                        act.rate.i = act.rate.i,
                        inf.prob.q = inf.prob.q,
                        act.rate.q = act.rate.q,
                        prog.rate = prog.rate,
                        quar.rate = quar.rate,
                        hosp.rate = hosp.rate,
                        disch.rate = disch.rate,
                        rec.rate = rec.rate,
                        arec.rate = arec.rate,
                        prog.dist.scale = prog.dist.scale,
                        prog.dist.shape = prog.dist.shape,
                        quar.dist.scale = quar.dist.scale,
                        quar.dist.shape = quar.dist.shape,
                        hosp.dist.scale = hosp.dist.scale,
                        hosp.dist.shape = hosp.dist.shape,
                        disch.dist.scale = disch.dist.scale,
                        disch.dist.shape = disch.dist.shape,
                        rec.dist.scale = rec.dist.scale,
                        rec.dist.shape = rec.dist.shape,
                        arec.dist.scale = arec.dist.scale,
                        arec.dist.shape = arec.dist.shape,
                        fat.rate.base = fat.rate.base,
                        hosp.cap = hosp.cap,
                        fat.rate.overcap = fat.rate.overcap,
                        fat.tcoeff = fat.tcoeff,
                        vital = vital,
                        a.rate = a.rate,
                        a.prop.e = a.prop.e,
                        a.prop.i = a.prop.i,
                        a.prop.q = a.prop.q,
                        ds.rate = ds.rate,
                        de.rate = de.rate,
                        di.rate = di.rate,
                        dq.rate = dq.rate,
                        dh.rate = dh.rate,
                        dr.rate = dr.rate)

    sim <- seiqhrf(param = param, init = init, control = control, seed = seed)
    sim_df <- as.data.frame(sim, out=out)
    class(sim) <- "icm"

    return(list(sim=sim, df=sim_df))
}
