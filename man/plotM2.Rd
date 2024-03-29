% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotM2}
\alias{plotM2}
\title{Alias for \code{plotPredMort()}}
\usage{
plotM2(
  object,
  species = NULL,
  time_range,
  all.sizes = FALSE,
  highlight = NULL,
  return_data = FALSE,
  ...
)
}
\arguments{
\item{object}{An object of class \linkS4class{MizerSim} or
\linkS4class{MizerParams}.}

\item{species}{The species to be selected. Optional. By default all target
species are selected. A vector of species names, or a
numeric vector with the species indices, or a logical vector indicating for
each species whether it is to be selected (TRUE) or not.}

\item{time_range}{The time range (either a vector of values, a vector of min
and max time, or a single value) to average the abundances over. Default is
the final time step. Ignored when called with a \linkS4class{MizerParams}
object.}

\item{all.sizes}{If TRUE, then predation mortality is plotted also for sizes
outside a species' size range. Default FALSE.}

\item{highlight}{Name or vector of names of the species to be highlighted.}

\item{return_data}{A boolean value that determines whether the formatted data
used for the plot is returned instead of the plot itself. Default value is FALSE}

\item{...}{Other arguments (currently unused)}
}
\value{
A ggplot2 object, unless \code{return_data = TRUE}, in which case a data
frame with the three variables 'w', 'value', 'Species' is returned.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
An alias provided for backward compatibility with mizer version <= 1.0
}
\examples{
\donttest{
params <-  NS_params
sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
plotPredMort(sim)
plotPredMort(sim, time_range = 10:20)

# Returning the data frame
fr <- plotPredMort(sim, return_data = TRUE)
str(fr)
}
}
\seealso{
\link{plotting_functions},  \code{\link[=getPredMort]{getPredMort()}}

Other plotting functions: 
\code{\link{animateSpectra}()},
\code{\link{plot,MizerSim,missing-method}},
\code{\link{plotBiomass}()},
\code{\link{plotDiet}()},
\code{\link{plotFMort}()},
\code{\link{plotFeedingLevel}()},
\code{\link{plotGrowthCurves}()},
\code{\link{plotSpectra}()},
\code{\link{plotYield}()},
\code{\link{plotYieldGear}()},
\code{\link{plotting_functions}}
}
\concept{deprecated}
