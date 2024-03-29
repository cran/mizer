% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotYield}
\alias{plotYield}
\alias{plotlyYield}
\title{Plot the total yield of species through time}
\usage{
plotYield(
  sim,
  sim2,
  species = NULL,
  total = FALSE,
  log = TRUE,
  highlight = NULL,
  return_data = FALSE,
  ...
)

plotlyYield(
  sim,
  sim2,
  species = NULL,
  total = FALSE,
  log = TRUE,
  highlight = NULL,
  ...
)
}
\arguments{
\item{sim}{An object of class \linkS4class{MizerSim}}

\item{sim2}{An optional second object of class \linkS4class{MizerSim}. If
this is provided its yields will be shown on the same plot in bolder lines.}

\item{species}{The species to be selected. Optional. By default all target
species are selected. A vector of species names, or a
numeric vector with the species indices, or a logical vector indicating for
each species whether it is to be selected (TRUE) or not.}

\item{total}{A boolean value that determines whether the total over all
species in the system is plotted as well. Note that even if the plot
only shows a selection of species, the total is including all species.
Default is FALSE.}

\item{log}{Boolean whether yield should be plotted on a logarithmic axis.
Defaults to true.}

\item{highlight}{Name or vector of names of the species to be highlighted.}

\item{return_data}{A boolean value that determines whether the formatted data
used for the plot is returned instead of the plot itself. Default value is FALSE}

\item{...}{Other arguments (currently unused)}
}
\value{
A ggplot2 object, unless \code{return_data = TRUE}, in which case a data
frame with the three variables 'Year', 'Yield', 'Species' is returned.
}
\description{
After running a projection, the total yield of each species across all
fishing gears can be plotted against time. The yield is obtained with
\code{\link[=getYield]{getYield()}}.
}
\examples{
\donttest{
params <- NS_params
sim <- project(params, effort = 1, t_max = 20, t_save = 0.2, progress_bar = FALSE)
plotYield(sim)
plotYield(sim, species = c("Cod", "Herring"), total = TRUE)

# Comparing with yield from twice the effort
sim2 <- project(params, effort=2, t_max=20, t_save = 0.2, progress_bar = FALSE)
plotYield(sim, sim2, species = c("Cod", "Herring"), log = FALSE)

# Returning the data frame
fr <- plotYield(sim, return_data = TRUE)
str(fr)
}
}
\seealso{
\link{plotting_functions},  \code{\link[=getYield]{getYield()}}

Other plotting functions: 
\code{\link{animateSpectra}()},
\code{\link{plot,MizerSim,missing-method}},
\code{\link{plotBiomass}()},
\code{\link{plotDiet}()},
\code{\link{plotFMort}()},
\code{\link{plotFeedingLevel}()},
\code{\link{plotGrowthCurves}()},
\code{\link{plotPredMort}()},
\code{\link{plotSpectra}()},
\code{\link{plotYieldGear}()},
\code{\link{plotting_functions}}
}
\concept{plotting functions}
