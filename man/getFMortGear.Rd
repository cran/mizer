% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rate_functions.R
\name{getFMortGear}
\alias{getFMortGear}
\title{Get the fishing mortality by time, gear, species and size}
\usage{
getFMortGear(object, effort, time_range)
}
\arguments{
\item{object}{A \code{MizerParams} object or a \code{MizerSim} object.}

\item{effort}{The effort for each fishing gear. See notes below.}

\item{time_range}{Subset the returned fishing mortalities by time. The time
range is either a vector of values, a vector of min and max time, or a
single value. Default is the whole time range. Only used if the
\code{object} argument is of type \code{MizerSim}.}
}
\value{
An array. If the effort argument has a time dimension, or a
\code{MizerSim} is passed in, the output array has four dimensions (time x
gear x species x size). If the effort argument does not have a time
dimension (i.e. it is a vector or a single numeric), the output array has
three dimensions (gear x species x size).
}
\description{
Calculates the fishing mortality rate \eqn{F_{g,i,w}} by gear, species and
size and possibly time (in units 1/year).
}
\note{
Here: fishing mortality = catchability x selectivity x effort.

The \code{effort} argument is only used if a \code{MizerParams} object is
passed in. The \code{effort} argument can be a two dimensional array (time x
gear), a vector of length equal to the number of gears (each gear has a
different effort that is constant in time), or a single numeric value (each
gear has the same effort that is constant in time). The order of gears in the
\code{effort} argument must be the same the same as in the \code{MizerParams}
object. If the \code{effort} argument is not supplied, its value is taken
from the \verb{@initial_effort} slot in the params object.

If the object argument is of class \code{MizerSim} then the effort slot of
the \code{MizerSim} object is used and the \code{effort} argument is not
used.
}
\examples{
\donttest{
params <-NS_params
# Get the fishing mortality when effort is constant
# for all gears and time:
getFMortGear(params, effort = 1)
# Get the fishing mortality when effort is different
# between the four gears but constant in time:
getFMortGear(params, effort = c(0.5, 1, 1.5, 0.75))
# Get the fishing mortality when effort is different
# between the four gears and changes with time:
effort <- array(NA, dim = c(20, 4))
effort[, 1] <- seq(from=0, to = 1, length = 20)
effort[, 2] <- seq(from=1, to = 0.5, length = 20)
effort[, 3] <- seq(from=1, to = 2, length = 20)
effort[, 4] <- seq(from=2, to = 1, length = 20)
getFMortGear(params, effort = effort)
# Get the fishing mortality using the effort already held in a MizerSim object.
sim <- project(params, t_max = 20, effort = 0.5)
getFMortGear(sim)
getFMortGear(sim, time_range = c(10, 20))
}

}
\seealso{
Other rate functions: 
\code{\link{getEGrowth}()},
\code{\link{getERepro}()},
\code{\link{getEReproAndGrowth}()},
\code{\link{getEncounter}()},
\code{\link{getFMort}()},
\code{\link{getFeedingLevel}()},
\code{\link{getMort}()},
\code{\link{getPredMort}()},
\code{\link{getPredRate}()},
\code{\link{getRDD}()},
\code{\link{getRDI}()},
\code{\link{getRates}()},
\code{\link{getResourceMort}()}
}
\concept{rate functions}
