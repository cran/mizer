% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_methods.R
\name{getMeanMaxWeight}
\alias{getMeanMaxWeight}
\title{Calculate the mean maximum weight of the community}
\usage{
getMeanMaxWeight(sim, species = NULL, measure = "both", ...)
}
\arguments{
\item{sim}{A \linkS4class{MizerSim} object}

\item{species}{The species to be selected. Optional. By default all target
species are selected. A vector of species names, or a
numeric vector with the species indices, or a logical vector indicating for
each species whether it is to be selected (TRUE) or not.}

\item{measure}{The measure to return. Can be 'numbers', 'biomass' or 'both'}

\item{...}{
  Arguments passed on to \code{\link[=get_size_range_array]{get_size_range_array}}
  \describe{
    \item{\code{min_w}}{Smallest weight in size range. Defaults to smallest weight in
the model.}
    \item{\code{max_w}}{Largest weight in size range. Defaults to largest weight in the
model.}
    \item{\code{min_l}}{Smallest length in size range. If supplied, this takes
precedence over \code{min_w}.}
    \item{\code{max_l}}{Largest length in size range. If supplied, this takes precedence
over \code{max_w}.}
  }}
}
\value{
Depends on the \code{measure} argument. If \code{measure = “both”}
  then you get a matrix with two columns, one with values by numbers, the
other with values by biomass at each saved time step. If \code{measure =
  “numbers”} or \code{“biomass”} you get a vector of the respective values at
each saved time step.
}
\description{
Calculates the mean maximum weight of the community through time. This can be
calculated by numbers or biomass. The calculation is the sum of the w_max *
abundance of each species, divided by the total abundance community, where
abundance is either in biomass or numbers. You can specify minimum and
maximum weight or length range for the species. Lengths take precedence over
weights (i.e. if both min_l and min_w are supplied, only min_l will be used).
You can also specify the species to be used in the calculation.
}
\examples{
mmw <- getMeanMaxWeight(NS_sim)
years <- c("1967", "2010")
mmw[years, ]
getMeanMaxWeight(NS_sim, species=c("Herring","Sprat","N.pout"))[years, ]
getMeanMaxWeight(NS_sim, min_w = 10, max_w = 5000)[years, ]
}
\seealso{
Other functions for calculating indicators: 
\code{\link{getCommunitySlope}()},
\code{\link{getMeanWeight}()},
\code{\link{getProportionOfLargeFish}()}
}
\concept{functions for calculating indicators}
\concept{summary_function}
