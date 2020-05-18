#' Define prior distributions to stochastic parameters
#'
#' @param param Model parameters, returned by \code{\link{param_seiqhrf}}. 
#' @param pars A vector of names of random parameters. 
#' @param prior.dist Probability prior distributions of the parameters in \code{par}. 
#' @param prior.pars Distribution parameters of the probability distributions in \code{prior.dist}. 
#' @param shift_only In effect only for parameters (in \code{param}) that have length longer than 1. 
#'    If TRUE, only parameter used at the first time point is simulated from prior distribution, 
#'    those used for the latter time points are shifted from the input values in param by the same amount as the first time point. 
#'   
#' 
#' @export

select_prior <- function(param, pars, prior.dist = NULL, prior.pars = list(), shift_only = FALSE){
    
    No_pars <- length(pars)
    ret <- vector("list", No_pars)
    names(ret) <- pars
    
    ### take care of empty dist and dist.par 
    if(is.null(prior.dist)){
       if(length(prior.pars) > 0) warning("prior.dist is not proveded, so prior.pars will be ignored")
       prior.dist <- vector("character", No_pars)  
       prior.pars <- vector("list", No_pars)
    }else{
        if(length(prior.pars) == 0) stop("If prior.dist is provided, prior.par needs to be specified as well")
    }
    
    ### Check inputs
    if(any(pars %in% names(param) == FALSE)) stop("All elements in pars need to specified in param")
    if(length(prior.dist) != No_pars) stop("pars and prior.dist must have the same length")
    if(length(prior.pars) != No_pars || !("list" %in% class(prior.pars))) stop("prior.pars need to be a list of the same length as pars")
    
    ### Check and organize into ret
    for(i in 1:No_pars){
        
        par <- pars[i]
        dist_tmp <- prior.dist[i]
        dist_par_tmp <- prior.pars[[i]]
        
        if(shift_only){
            mu <- param[[par]][1]
        }else{
            mu <- param[[par]]
        }
        
        
        par_length <- length(mu)
        
        ### check user specified prior.dist and dist.par
        if(dist_tmp != ""){
            
            if(!(dist_tmp %in% c("gaussian", "beta"))) stop("prior.dist needs to be either gaussian or beta")
            if(grepl("prob", par) && dist_tmp != "beta") stop(paste0("Parameter ", par, " seems to be a probability, please assign prior.dist = NULL or prior.dist = 'beta'"))
            if(!("list" %in% class(dist_par_tmp))) stop(paste0("The ", i, "th entry in dist.par is not a list"))
            
            if(dist_tmp == "gaussian"){
                if(is.null(dist_par_tmp[["sd"]])) stop("dist.par corresponding to prior.dist == gaussian requires at least 'sd' value")
                if(is.null(dist_par_tmp[["mean"]])) dist_par_tmp[["mean"]] <- mu
                if(length(dist_par_tmp[["sd"]]) == 1) dist_par_tmp[["sd"]] <- rep(dist_par_tmp[["sd"]], par_length)
            }
            if(dist_tmp == "beta"){
                if(is.null(dist_par_tmp[["a"]]) && is.null(dist_par_tmp[["b"]])) stop("If prior.dist == beta, at least one of 'a' and 'b' needs to be specified")
                if(is.null(dist_par_tmp[["a"]])){
                    if(length(dist_par_tmp[["b"]]) == 1) dist_par_tmp[["b"]] <- rep(dist_par_tmp[["b"]], par_length)
                    dist_par_tmp[["a"]] <- dist_par_tmp[["b"]]*mu/(1 - mu)
                }else{
                    if(length(dist_par_tmp[["a"]]) == 1) dist_par_tmp[["a"]] <- rep(dist_par_tmp[["a"]], par_length)
                    dist_par_tmp[["b"]] <- dist_par_tmp[["a"]]*(1 - mu)/mu
                }
            }
            
                
       ### Assign default prior.dist and dist.par   
        }else{
            if(!is.null(dist_par_tmp)) warning(paste0("The prior.dist of ", par, " is not specified, so dist.par will be ignored"))
            if(grepl("act", par)){
                dist_tmp <- "gaussian"
                dist_par_tmp[["mean"]] <- mu
                dist_par_tmp[["sd"]] <- rep(0.1, par_length)
            }else{
                dist_tmp <- "beta"
                dist_par_tmp[["a"]] <- rep(50, par_length)
                dist_par_tmp[["b"]] <- 50*(1 - mu)/mu
            }  
        }
       
       
       ### save as ret 
        ret[[par]] <- list("dist" = dist_tmp, "par" = dist_par_tmp)
    }
    
    ### return
    ret$shift_only <- shift_only
    class(ret) <- "prior"
    ret
}




#' Print important information from specified prior.
#' @param x A prior object returned by \code{\link{select_prior}}.
#' @param ... Additional parameters.
#' @export
print.prior <- function(x, ...){
    
    print_pars <- function(par1, par2){
        print1 <- NULL
        print2 <- NULL
        if(length(par1) > 1){
            for(i in 1:min(3, length(par1))){
                print1 <- paste0(print1, round(par1[i], 3), ", ")
                print2 <- paste0(print2, round(par2[i], 3), ", ")
            }
            if(length(par1) > 3){
                print1 <- paste0(print1, "...")
                print2 <- paste0(print2, "...")
            }
        }else return(c(round(par1, 3), round(par2, 3)))
        return(c(print1, print2))
    }
    
    cat("Prior distributions for SEIQHRF Parameters")
    cat("\n===========================================\n")
    
    for(par in names(x)[names(x) != "shift_only"]){
        if(x[[par]]$dist == "beta"){
            a <- x[[par]]$par[["a"]]
            b <- x[[par]]$par[["b"]]
            prints <- print_pars(a, b)
            prints_inference <- print_pars(a/(a + b), sqrt(a*b/((a+b)^2*(a + b + 1))))
            cat(par, ": Beta( a =", prints[1], ", b =", prints[2], "), mean =", prints_inference[1], ", s.e =", prints_inference[2], "\n")
        } 
        if(x[[par]]$dist == "gaussian"){
            prints <- print_pars(x[[par]]$par[["mean"]], x[[par]]$par[["sd"]])
            cat(par, ": N( mean =", prints[1], ", sd =", prints[2], ")\n")  
        } 
    }
    cat("===========================================\n")
}

# Internal 
#' @import ggplot2
extract_dist <- function(par, par_name, at = NULL){
    
    dist <- par$dist
    No_pars <- length(par$par[[1]])
    
    if(No_pars > 10 && is.null(at)){
        at <- round(seq(1, No_pars, No_pars/10))
        No_pars <- length(at)
    } else{
        at <- c(1:No_pars)
    }
    
    if(dist == "beta"){
        par_names <- c("a", "b")
        dfun <- stats::dbeta
        qfun <- stats::qbeta
        par1 <- par$par[["a"]][at]
        par2 <- par$par[["b"]][at]
    } 
    if(dist == "gaussian"){
        par_names <- c("mean", "sd")
        dfun <- stats::dnorm 
        qfun <- stats::qnorm
        par1 <- par$par[["mean"]][at]
        par2 <- par$par[["sd"]][at]
    } 
    
    low <- min(apply(rbind(par1, par2), 2, function(x, fun) fun(0.01, x[1], x[2]), fun = qfun))
    up <- max(apply(rbind(par1, par2), 2, function(x, fun) fun(0.99, x[1], x[2]), fun = qfun))
    
    x <- seq(low, up, (up - low)/500)
    plot_df <- NULL
    for(i in 1:No_pars){
        plot_df <- rbind(plot_df, data.frame(x = x, density = dfun(x, par1[i], par2[i]), Step = as.factor(at[i]))) 
    }
    
    ggplot(plot_df, aes(x = x, y = density, color = Step)) + 
        geom_line(size = 0.8, alpha = 0.8) + ggtitle(par_name)
    
}



#' Plot theoretical probability density of specified priors.
#' @param x A prior object returned by \code{\link{select_prior}}.
#' @param at A vector of time points to plot when an element in x has length greater than 1. 
#' @param ... Additional parameters.
#' 
#' @importFrom stats dbeta
#' @importFrom stats qbeta
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom gridExtra grid.arrange
#' @export
plot.prior <- function(x, at = NULL, ...){
    
    ret <- list()
    for(i in 1:(length(x) - 1)){
        ret[[i]] <- extract_dist(x[[i]], names(x)[i], at)
    }
    ret 
}