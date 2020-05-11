#' Generating time-dependent parameters
#'
#' Function to generate parameter values for the time range of the simulation
#' that can change over time.
#'
#' @param nstep Number of time steps to generate parameter values for.
#' @param vals List of parameter values to include over the nsteps.
#' @param timing List of the step numbers at which to start changes, with the
#'        last number reflecting when to hit the last parameter value in vals.
#'
#' @return list of parameter values for length t
#'
#' @importFrom utils tail
#' @export
vary_param <- function(nstep, vals, timing) {

    stopifnot(length(vals) == length(timing))
    y <- list()

    for(t in seq(1:nstep)){
        if(t <= timing[1]){              # If before first jump set to val[1]
            y.t <- vals[1]
        }else if(t > utils::tail(timing, n=1)){ # If after last jump set to val[-1]
            y.t <- utils::tail(vals, n=1)
        }else{                           # If intermediate step, calculate...
            for(j in (2:length(timing))){
                if(t > timing[j - 1] & t <= timing[j]){
                    start <- vals[j - 1]
                    end <- vals[j]
                    start_t <- timing[j - 1]
                    end_t <- timing[j]
                    y.t <- start - (t-start_t)*(start - end)/(end_t - start_t)
                }
            }
        }
        y <- append(y, y.t)
    }

    return(unlist(y))

}
