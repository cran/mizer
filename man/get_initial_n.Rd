% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_initial_n.R
\name{get_initial_n}
\alias{get_initial_n}
\title{Calculate initial population abundances}
\usage{
get_initial_n(params, n0_mult = NULL, a = 0.35)
}
\arguments{
\item{params}{The model parameters. An object of type
\linkS4class{MizerParams}.}

\item{n0_mult}{Multiplier for the abundance at size 0. Default value is
kappa/1000.}

\item{a}{A parameter with a default value of 0.35.}
}
\value{
A matrix (species x size) of population abundances.
}
\description{
This function uses the model parameters and other parameters to calculate
initial values for the species number densities. These initial
abundances are currently quite arbitrary and not close to the steady state.
We intend to improve this in the future.
}
\examples{
init_n <- get_initial_n(NS_params)
}
\concept{helper}
