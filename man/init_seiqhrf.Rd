% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PREP_init.R
\name{init_seiqhrf}
\alias{init_seiqhrf}
\title{Initial Conditions for Stochastic Individual Contact Models}
\usage{
init_seiqhrf(
  s.num = 9997,
  e.num = 0,
  i.num = 3,
  q.num = 0,
  h.num = 0,
  r.num = 0,
  f.num = 0,
  ...
)
}
\arguments{
\item{s.num}{Initial number of *S compartment individuals in
the simulated population. An overall population of 10,000 is a good
compromise. A set of models will still take several minutes or more
to run, in parallel.}

\item{e.num}{Initial number of E compartment individuals in
the simulated population.}

\item{i.num}{Initial number of I compartment individuals in
the simulated population.}

\item{q.num}{Initial number of Q compartment individuals in
the simulated population.}

\item{h.num}{Initial number of H compartment individuals in
the simulated population.}

\item{r.num}{Initial number of R compartment individuals in
the simulated population.}

\item{f.num}{Initial number of F compartment individuals in
the simulated population.}

\item{...}{Additional initial conditions passed to model.}
}
\description{
Sets the initial conditions for stochastic individual contact
             models simulated with \code{icm}.
}
\details{
The initial conditions for a model solved with \code{\link{icm}} should be
input into the \code{init.icm} function. This function handles initial
conditions for both base models and original models using new modules. For
an overview of initial conditions for base ICM class models, consult the
\href{http://statnet.github.io/tut/BasicICMs.html}{Basic ICMs} tutorial.
}
\seealso{
Use \code{\link{param.icm}} to specify model parameters and
         \code{\link{control.icm}} to specify the control settings. Run the
         parameterized model with \code{\link{icm}}.
}
\keyword{parameterization}
