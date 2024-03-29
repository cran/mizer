% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_methods.R
\name{getProportionOfLargeFish}
\alias{getProportionOfLargeFish}
\title{Calculate the proportion of large fish}
\usage{
getProportionOfLargeFish(
  sim,
  species = NULL,
  threshold_w = 100,
  threshold_l = NULL,
  biomass_proportion = TRUE,
  ...
)
}
\arguments{
\item{sim}{A \linkS4class{MizerSim} object}

\item{species}{The species to be selected. Optional. By default all target
species are selected. A vector of species names, or a
numeric vector with the species indices, or a logical vector indicating for
each species whether it is to be selected (TRUE) or not.}

\item{threshold_w}{the size used as the cutoff between large and small fish.
Default value is 100.}

\item{threshold_l}{the size used as the cutoff between large and small fish.}

\item{biomass_proportion}{a boolean value. If TRUE the proportion calculated
is based on biomass, if FALSE it is based on numbers of individuals.
Default is TRUE.}

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
A vector containing the proportion of large fish through time
}
\description{
Calculates the proportion of large fish through time in the \code{MizerSim}
class within user defined size limits. The default option is to use the whole
size range. You can specify minimum and maximum size ranges for the species
and also the threshold size for large fish. Sizes can be expressed as weight
or size. Lengths take precedence over weights (i.e. if both min_l and min_w
are supplied, only min_l will be used). You can also specify the species to
be used in the calculation. This function can be used to calculate the Large
Fish Index. The proportion is based on either abundance or biomass.
}
\examples{
lfi <- getProportionOfLargeFish(NS_sim, min_w = 10, max_w = 5000, 
                                threshold_w = 500)
years <- c("1972", "2010")
lfi[years]
getProportionOfLargeFish(NS_sim)[years]
getProportionOfLargeFish(NS_sim, species=c("Herring","Sprat","N.pout"))[years]
getProportionOfLargeFish(NS_sim, min_w = 10, max_w = 5000)[years]
getProportionOfLargeFish(NS_sim, min_w = 10, max_w = 5000,
    threshold_w = 500, biomass_proportion = FALSE)[years]
}
\seealso{
Other functions for calculating indicators: 
\code{\link{getCommunitySlope}()},
\code{\link{getMeanMaxWeight}()},
\code{\link{getMeanWeight}()}
}
\concept{functions for calculating indicators}
\concept{summary_function}
