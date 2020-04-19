#' Sample value from Weibull distribution
#'
#' Function to sample values from a Weibull distribution with parameters shape
#' and rate
#'
#' @param vecTimeSinceExp ?
#' @param scale scale parameter value
#' @param shape shape parameter value
#' 
#' @return Value of coefficient of variation for vector
#' @export
cum_discr_si <- function(vecTimeSinceExp, scale, shape) {
    vlen <- length(vecTimeSinceExp)
    if (vlen > 0) {
        probVec <- numeric(vlen)
        for (p in 1:vlen) {
            probVec[p] <- pweibull(vecTimeSinceExp[p], shape=shape, scale=scale)
        }
    } else {
        probVec <- 0
    }
    return(probVec)
}
