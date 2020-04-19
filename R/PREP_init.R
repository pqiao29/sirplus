#' @title Initial Conditions for Stochastic Individual Contact Models
#'
#' @description Sets the initial conditions for stochastic individual contact
#'              models simulated with \code{icm}.
#'
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
#' @param ... Additional initial conditions passed to model.
#'
#' @details
#' The initial conditions for a model solved with \code{\link{icm}} should be
#' input into the \code{init.icm} function. This function handles initial
#' conditions for both base models and original models using new modules. For
#' an overview of initial conditions for base ICM class models, consult the
#' \href{http://statnet.github.io/tut/BasicICMs.html}{Basic ICMs} tutorial.
#'
#' @seealso Use \code{\link{param.icm}} to specify model parameters and
#'          \code{\link{control.icm}} to specify the control settings. Run the
#'          parameterized model with \code{\link{icm}}.
#'
#' @keywords parameterization
#'
#' @export
#'

init_seiqhrf <- function(s.num = 9997, e.num=0, i.num = 3, q.num=0, h.num=0, r.num = 0, f.num = 0, ...) {

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
    
    
    ## Output
    class(p) <- c("init.seiqhrf", "list")
    return(p)
}