% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steady.R
\name{steady}
\alias{steady}
\title{Set initial values to a steady state for the model}
\usage{
steady(
  params,
  t_max = 100,
  t_per = 1.5,
  dt = 0.1,
  tol = 0.1 * dt,
  return_sim = FALSE,
  preserve = c("reproduction_level", "erepro", "R_max"),
  progress_bar = TRUE
)
}
\arguments{
\item{params}{A \linkS4class{MizerParams} object}

\item{t_max}{The maximum number of years to run the simulation. Default is 100.}

\item{t_per}{The simulation is broken up into shorter runs of \code{t_per} years,
after each of which we check for convergence. Default value is 1.5. This
should be chosen as an odd multiple of the timestep \code{dt} in order to be
able to detect period 2 cycles.}

\item{dt}{The time step to use in \code{project()}.}

\item{tol}{The simulation stops when the relative change in the egg
production RDI over \code{t_per} years is less than \code{tol} for every species.}

\item{return_sim}{If TRUE, the function returns the MizerSim object holding
the result of the simulation run, saved at intervals of \code{t_per}. If FALSE (default) the function returns
a MizerParams object with the "initial" slots set to the steady state.}

\item{preserve}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
Specifies whether the \code{reproduction_level} should be preserved (default)
or the maximum reproduction rate \code{R_max} or the reproductive
efficiency \code{erepro}. See \code{\link[=setBevertonHolt]{setBevertonHolt()}} for an explanation
of the \code{reproduction_level}.}

\item{progress_bar}{A shiny progress object to implement a progress bar in a
shiny app. Default FALSE.}
}
\value{
A MizerParams or a MizerSim object
}
\description{
The steady state is found by running the dynamics while keeping reproduction,
resource and other components constant until the size spectra no longer
change much (or until time \code{t_max} is reached, if earlier).
}
\details{
If the model use Beverton-Holt reproduction then the reproduction parameters
are set to values that give the level of reproduction observed in that
steady state. The \code{preserve} argument can be used to specify which of the
reproduction parameters should be preserved.
}
\examples{
\donttest{
params <- newTraitParams()
species_params(params)$gamma[5] <- 3000
params <- steady(params)
plotSpectra(params)
}
}
