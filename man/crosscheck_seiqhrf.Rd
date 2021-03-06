% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PREP_crosscheck.R
\name{crosscheck_seiqhrf}
\alias{crosscheck_seiqhrf}
\title{Cross Checking of Inputs for Stochastic Individual Contact Models}
\usage{
crosscheck_seiqhrf(param, init, control)
}
\arguments{
\item{param}{An \code{EpiModel} object of class \code{param.seiqhrf}.}

\item{init}{An \code{EpiModel} object of class \code{init.icm}.}

\item{control}{A list returned by \code{\link{control_seiqhrf}}.}
}
\value{
This function returns no objects.
}
\description{
This function checks that the three parameter lists from
             \code{\link{param_seiqhrf}}, \code{\link{init.icm}}, and
             \code{\link{control_seiqhrf}} are consistent.
}
\keyword{internal}
