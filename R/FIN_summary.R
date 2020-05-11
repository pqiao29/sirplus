#' summary function
#'
#' Function to extract mean and standard deviation of all compartments 
#' for SEIQHRF model
#'
#' @param object seiqhrf object, returned by \code{\link{seiqhrf}}.
#' @param ... Additional parameters.
#'
#' @return A list of compartments, each compartment is a list of three: 
#' \itemize{
#' \item \code{mean}: (t-dimensional vector) Mean of simulated counts of the compartment;
#' \item \code{sd}: (t-dimensional vector) Standard deviation of simulated counts of the compartment;
#' \item \code{CI}: (t by 2 matrix) 95\% confidence intervals of simulated counts of the compartment (assumed Gaussian),
#' \item \code{qntCI}: (t by 2 matrix) 95\% quantile confidence intervals,
#' }
#' where t is the total number of time points in the simulation (i.e. object$nstep).
#'
#' @importFrom stats setNames
#' @export
summary.seiqhrf <- function(object, ...){
    
    nsteps <- object$control$nsteps
    c_names <- c("s.num", "e.num", "i.num", "q.num", "h.num", "r.num", "f.num")
    
    #### extract mean and sd
    df_mean <- as.data.frame(object, out = "mean")
    sub_cols <- names(df_mean)[ order(match(names(df_mean), c_names))][1:length(c_names)]
    df_mean <- df_mean[, sub_cols]  
    
    df_sd <- as.data.frame(object, out = "sd")[, sub_cols]

    df_qnt_left <- as.data.frame(object, out = "qnt", qval = .025)[, sub_cols]
    df_qnt_right <- as.data.frame(object, out = "qnt", qval = .975)[, sub_cols]
    
    #### prepare result 
    res <- stats::setNames(vector("list", length(c_names)), c_names)

    for(prev_no in c_names){

        prev_mean <- df_mean[[prev_no]]
        prev_sd <- df_sd[[prev_no]]
        prev_ci <- cbind(prev_mean - 1.96*prev_sd, prev_mean + 1.96*prev_sd)
        prev_qci <- cbind(df_qnt_left[[prev_no]], df_qnt_right[[prev_no]])
        names(prev_mean) <- c(1:nsteps)
        names(prev_sd) <- c(1:nsteps)
        rownames(prev_ci) <- c(1:nsteps)
        rownames(prev_qci) <- c(1:nsteps)
        
        res[[prev_no]]$"mean" <- prev_mean
        res[[prev_no]]$"sd" <- prev_sd
        res[[prev_no]]$"CI" <- prev_ci
        res[[prev_no]]$"qntCI" <- prev_qci
        
    }
    
    class(res) <- "summary.seiqhrf"
    res
}

