% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibrate.R
\name{scaleModel}
\alias{scaleModel}
\title{Change scale of the model}
\usage{
scaleModel(params, factor)
}
\arguments{
\item{params}{A MizerParams object}

\item{factor}{The factor by which the scale is multiplied}
}
\value{
The rescaled MizerParams object
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

The abundances in mizer and some rates depend on the size of the area to
which they refer. So they could be given per square meter or per square
kilometer or for an entire study area or any other choice of yours. This
function allows you to change the scale of the model by automatically
changing the abundances and rates accordingly.
}
\details{
If you rescale the model by a factor \eqn{c} then this function makes the
following rescalings in the params object:
\itemize{
\item The initial abundances are rescaled by \eqn{c}.
\item The search volume is rescaled by \eqn{1/c}.
\item The resource carrying capacity is rescaled by \eqn{c}
\item The maximum reproduction rate \eqn{R_{max}} is rescaled by
\eqn{c}.
}
The effect of this is that the dynamics of the rescaled model are identical
to those of the unscaled model, in the sense that it does not matter whether
one first calls \code{\link[=scaleModel]{scaleModel()}} and then runs a simulation with
\code{\link[=project]{project()}} or whether one first runs a simulation and then rescales the
resulting abundances.

Note that if you use non-standard resource dynamics or other components then
you may need to rescale additional parameters that appear in those dynamics.

In practice you will need to use some observations to set the scale for your
model. If you have biomass observations you can use \code{\link[=calibrateBiomass]{calibrateBiomass()}},
if you have yearly yields you can use \code{\link[=calibrateYield]{calibrateYield()}}.
}
