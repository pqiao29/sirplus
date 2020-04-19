#' @title Cross Checking of Inputs for Stochastic Individual Contact Models
#'
#' @description This function checks that the three parameter lists from
#'              \code{\link{param_seiqhrf}}, \code{\link{init.icm}}, and
#'              \code{\link{control_seiqhrf}} are consistent.
#'
#' @param param An \code{EpiModel} object of class \code{param.seiqhrf}.
#' @param init An \code{EpiModel} object of class \code{init.icm}.
#' @param control A list returned by \code{\link{control_seiqhrf}}.
#'
#' @return
#' This function returns no objects.
#'
#' @keywords internal
#'
crosscheck_seiqhrf <- function(param, init, control) {
    
    ### Main class check
    if (!inherits(param, "param.seiqhrf"))  stop("param must an object of class param.seiqhrf", call. = FALSE)
    if (!inherits(control, "control.seiqhrf"))  stop("control must an object of class control.seiqhrf", call. = FALSE)
        
        
    
    if (control$skip.check == FALSE) {
        
        ## Check that rec.rate is supplied for SIR models
        if (control$type %in% c("SIR", "SIS") && is.null(param$rec.rate)) {
                stop("Specify rec.rate in param.seiqhrf", call. = FALSE)
        }
        
        
        ## Check that paramets and init are supplied for SIR models
        if (control$type == "SIR" && is.null(init$r.num)) {
            stop("Specify r.num in init.icm", call. = FALSE)
        }
        
        
        ## Deprecated parameters
        bim <- grep(".FUN", names(control), value = TRUE)
        um <- which(grepl(".FUN", names(control)) & !(names(control) %in% bim))
        if (length(um) == 0 && !is.null(control$type)) {
            if (!is.null(param$trans.rate)) {
                stop("The trans.rate parameter is deprecated. Use the inf.prob ",
                     "parameter instead.", call. = FALSE)
            }
        }
        
        ## Check param (modified from infection.FUN and arrivals.FUN) -------------------------------------
        check_list <- c("act.rate.i", "inf.prob.i", "a.rate", "a.prop.e", "a.prop.i", "a.prop.q")
        if(control$type %in% c("SEIQHR", "SEIQHRF")) check_list<- c(check_list, "quar.rate", "disch.rate", "act.rate.e", "inf.prob.e", "act.rate.q", "inf.prob.q")
        check_idx <- which(names(param) %in% check_list)
        
        for(i in check_idx){
            param_length = length(param[[i]])
            if (param_length != 1 && param_length != control$nsteps) stop(paste0("Length of", names(param)[i], "must be 1 or the value of nsteps"))
        }
    }
    
    ## In-place assignment to update param and control
    on.exit(assign("param", param, pos = parent.frame()))
    on.exit(assign("control", control, pos = parent.frame()), add = TRUE)
}