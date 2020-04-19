#' @export
print.init.seiqhrf  <- function(x, ...) {
    
    cat("SEIQHRF Initial Conditions")
    cat("\n===========================\n")
    cat("\nUser specified control parameters:")
    cat("\n---------------------------\n")
    for(i in x$usr.specified){
        if (class(x[[i]]) %in% c("integer", "numeric") && length(x[[i]]) > 10) {
            cat(names(x)[i], "=", x[[i]][1:5], "...", fill = 80)
        } else if (class(x[[i]]) == "data.frame") {
            cat(names(x)[i], "= <data.frame>\n")
        } else if (class(x[[i]]) == "list") {
            cat(i, "= <list>\n")
        } else {
            cat(i, "=", x[[i]], fill = 80)
        }
    }
    cat("\nDefault control parameters:")
    cat("\n---------------------------\n")
    
    pToPrint <- names(x)[seq_along(x)]
    pToPrint<- setdiff(pToPrint, c(x$usr.specified, "usr.specified"))
    for (i in pToPrint) {
        if (class(x[[i]]) %in% c("integer", "numeric") && length(x[[i]]) > 10) {
            cat(names(x)[i], "=", x[[i]][1:5], "...", fill = 80)
        } else if (class(x[[i]]) == "data.frame") {
            cat(names(x)[i], "= <data.frame>\n")
        } else if (class(x[[i]]) == "list") {
            cat(i, "= <list>\n")
        } else {
            cat(i, "=", x[[i]], fill = 80)
        }
    }
    
    invisible()
}


#' @export
print.param.seiqhrf <- function(x, ...) {
    
    cat("SEIQHRF Parameters")
    cat("\n===========================\n")
    
    cat("\nUser specified control parameters:")
    cat("\n---------------------------\n")
    for(i in x$usr.specified){
        if (class(x[[i]]) %in% c("integer", "numeric") && length(x[[i]]) > 10) {
            cat(i, "=", x[[i]][1:5], "...", fill = 80)
        } else if (class(x[[i]]) == "data.frame") {
            cat(i, "= <data.frame>\n")
        } else if (class(x[[i]]) == "list") {
            cat(i, "= <list>\n")
        } else {
            cat(i, "=", x[[i]], fill = 80)
        }
    }
    
    cat("\nDefault control parameters:")
    cat("\n---------------------------\n")
    pToPrint <- names(x)[which(!(names(x) %in% c("vital")))]
    pToPrint<- setdiff(pToPrint, c(x$usr.specified, "usr.specified"))
    for (i in pToPrint) {
        if (class(x[[i]]) %in% c("integer", "numeric") && length(x[[i]]) > 10) {
            cat(i, "=", x[[i]][1:5], "...", fill = 80)
        } else if (class(x[[i]]) == "data.frame") {
            cat(i, "= <data.frame>\n")
        } else if (class(x[[i]]) == "list") {
            cat(i, "= <list>\n")
        } else {
            cat(i, "=", x[[i]], fill = 80)
        }
    }
    
    invisible()
}


#' @export
print.control.seiqhrf <- function(x, ...) {
    
    cat("SEIQHRF Control Settings")
    cat("\n===========================\n")
    cat("\nUser specified control parameters:")
    cat("\n---------------------------\n")
    for(i in x$usr.specified){
        cat(i, "=", x[[i]], fill = 80)
    }
    cat("\nDefault control parameters:")
    cat("\n---------------------------\n")
    pToPrint <- names(x)[which(!grepl(".FUN", names(x)))]
    pToPrint<- setdiff(pToPrint, c(x$usr.specified, "usr.specified"))
    for (i in pToPrint) {
        cat(i, "=", x[[i]], fill = 80)
    }
    cat("Base Modules:", x$bi.mods, fill = 80)
    if (length(x$user.mods) > 0) {
        cat("Extension Modules:", x$user.mods, fill = 80)
    }
    
    invisible()
}

#' @export
print.seiqhrf <- function(x, ...) {
    
    cat("SEIQHRF Model Simulation")
    cat("\n=======================")
    cat("\nModel class:", class(x))
    
    cat("\n\nSimulation Summary")
    cat("\n-----------------------")
    cat("\nModel type:", x$control$type)
    cat("\nNumber of simulations:", x$control$nsims)
    cat("\nNumber of time steps:", x$control$nsteps)
    
    cat("\n \nModel Output (variable names)")
    cat("\n-----------------------\n")
    cat(names(x$epi), fill = 60)
    
    invisible()
}

#' @export
print.summary.seiqhrf <- function(x, ...){
    
    cat("Peaks of SEIQHRF Model Simulation: ")
    cat("\n===============================\n")
    
    No_comps <- length(x)
    ret <- matrix(NA, No_comps, 2)
    for(i in 1:nrow(ret)){
        tmp_mean <- x[[i]]$mean
        tmp_mean[is.na(tmp_mean)] <- 0
        ret[i, 1] <- max(tmp_mean)
        ret[i, 2] <- which(x[[i]]$mean == ret[i, 1])[1]
    }
    colnames(ret) <- c("Max", "Time")
    rownames(ret) <- names(x)

   ret[, 1] <- format(ret[, 1], nsmall = 2, digits = 2)
   print(ret, quote = FALSE, right = TRUE)
}
