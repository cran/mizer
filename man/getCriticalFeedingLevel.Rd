% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rate_functions.R
\name{getCriticalFeedingLevel}
\alias{getCriticalFeedingLevel}
\title{Get critical feeding level}
\usage{
getCriticalFeedingLevel(params)
}
\arguments{
\item{params}{A MizerParams object}
}
\value{
A matrix (species x size) with the critical feeding level
}
\description{
The critical feeding level is the feeding level at which the food intake is
just high enough to cover the metabolic costs, with nothing left over for
growth or reproduction.
}
