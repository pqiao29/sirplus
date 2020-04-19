#' @title Extract Model Data for Stochastic Models
#'
#' @description This function extracts model simulations for objects of classes
#'              \code{seiqhrf} into a data frame using
#'              the generic \code{as.data.frame} function.
#'
#' @param x An \code{EpiModel} object of class \code{icm} or \code{netsim}.
#' @param out Data output to data frame: \code{"mean"} for row means across
#'        simulations, \code{"sd"} for row standard deviations across simulations,
#'        \code{"qnt"} for row quantiles at the level specified in \code{qval},
#'        or \code{"vals"} for values from individual simulations.
#' @param sim If \code{out="vals"}, the simulation number to output. If not
#'        specified, then data from all simulations will be output.
#' @param qval Quantile value required when \code{out="qnt"}.
#' @param row.names See \code{\link{as.data.frame.default}}.
#' @param optional See \code{\link{as.data.frame.default}}.
#' @param ...  See \code{\link{as.data.frame.default}}.
#'
#' @details
#' These methods work for both \code{seiqhrf} class models. The
#' available output includes time-specific means, standard deviations, quantiles,
#' and simulation values (compartment and flow sizes) from these stochastic model
#' classes. Means, standard deviations, and quantiles are calculated by taking the
#' row summary (i.e., each row of data is corresponds to a time step) across all
#' simulations in the model output.
#'
#' @method as.data.frame seiqhrf
#' @keywords extract
#' @export
#'
#' @examples
#' ## Stochastic SEIQHRF model
#' param <- param_seiqhrf()
#' init <- init_seiqhrf(s.num = 500, i.num = 1)
#' control <- control_seiqhrf(nsteps = 10, nsims = 3)
#' mod <- seiqhrf(init, control, param)
#'
#' # Default output all simulation runs, default to all in stacked data.frame
#' as.data.frame(mod)
#' as.data.frame(mod, sim = 2)
#'
#' # Time-specific means across simulations
#' as.data.frame(mod, out = "mean")
#'
#' # Time-specific standard deviations across simulations
#' as.data.frame(mod, out = "sd")
#'
#' # Time-specific quantile values across simulations
#' as.data.frame(mod, out = "qnt", qval = 0.25)
#' as.data.frame(mod, out = "qnt", qval = 0.75)
#'

as.data.frame.seiqhrf <- function(x, row.names = NULL, optional = FALSE,
                              out = "vals", sim, qval, ...) {
    
    df <- data.frame(time = 1:x$control$nsteps)
    nsims <- x$control$nsims
    
    if (out == "vals") {
        
        if (missing(sim)) {
            sim <- 1:nsims
        }
        if (max(sim) > nsims) {
            stop("Maximum sim is ", nsims, call. = FALSE)
        }
        
        for (j in sim) {
            if (j == min(sim)) {
                for (i in seq_along(x$epi)) {
                    if (nsims == 1) {
                        df[, i + 1] <- x$epi[[i]]
                    } else {
                        df[, i + 1] <- x$epi[[i]][, j]
                    }
                }
                df$sim <- j
            } else {
                tdf <- data.frame(time = 1:x$control$nsteps)
                for (i in seq_along(x$epi)) {
                    if (nsims == 1) {
                        tdf[, i + 1] <- x$epi[[i]]
                    } else {
                        tdf[, i + 1] <- x$epi[[i]][, j]
                    }
                }
                tdf$sim <- j
                df <- rbind(df, tdf)
            }
        }
        df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
        
    }
    
    ## Output means
    if (out == "mean") {
        if (nsims == 1) {
            for (i in seq_along(x$epi)) {
                df[, i + 1] <- x$epi[[i]]
            }
        }
        if (nsims > 1) {
            for (i in seq_along(x$epi)) {
                df[, i + 1] <- rowMeans(x$epi[[i]], na.rm = TRUE)
            }
        }
    }
    
    ## Output standard deviations
    if (out == "sd") {
        if (nsims == 1) {
            for (i in seq_along(x$epi)) {
                df[, i + 1] <- 0
            }
        }
        if (nsims > 1) {
            for (i in seq_along(x$epi)) {
                df[, i + 1] <- apply(x$epi[[i]], 1, sd, na.rm = TRUE)
            }
        }
    }
    
    if (out == "qnt") {
        if (missing(qval) || length(qval) > 1 || (qval > 1 | qval < 0)) {
            stop("Must specify qval as single value between 0 and 1", call. = FALSE)
        }
        for (i in seq_along(x$epi)) {
            df[, i + 1] <- apply(x$epi[[i]], 1, quantile, probs = qval,
                                 na.rm = TRUE, names = FALSE)
        }
    }
    
    if (out == "vals") {
        names(df)[3:ncol(df)] <- names(x$epi)
    } else {
        names(df)[2:ncol(df)] <- names(x$epi)
    }
    
    return(df)
}
