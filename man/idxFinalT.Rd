% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MizerSim-class.R
\name{idxFinalT}
\alias{idxFinalT}
\title{Time index at end of simulation}
\usage{
idxFinalT(sim)
}
\arguments{
\item{sim}{A MizerSim object}
}
\value{
An integer giving the index for extracting the
results for the final time step
}
\description{
Time index at end of simulation
}
\examples{
idx <- idxFinalT(NS_sim)
idx
# This coincides with
length(getTimes(NS_sim))
# and corresponds to the final time
getTimes(NS_sim)[idx]
# We can use this index to extract the result at the final time
identical(N(NS_sim)[idx, , ], finalN(NS_sim))
identical(NResource(NS_sim)[idx, ], finalNResource(NS_sim))
}
