% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resource_semichemostat.R
\name{resource_semichemostat}
\alias{resource_semichemostat}
\alias{balance_resource_semichemostat}
\title{Project resource using semichemostat model}
\usage{
resource_semichemostat(
  params,
  n,
  n_pp,
  n_other,
  rates,
  t,
  dt,
  resource_rate,
  resource_capacity,
  ...
)

balance_resource_semichemostat(params, resource_rate, resource_capacity)
}
\arguments{
\item{params}{A \link{MizerParams} object}

\item{n}{A matrix of species abundances (species x size)}

\item{n_pp}{A vector of the resource abundance by size}

\item{n_other}{A list with the abundances of other components}

\item{rates}{A list of rates as returned by \code{\link[=mizerRates]{mizerRates()}}}

\item{t}{The current time}

\item{dt}{Time step}

\item{resource_rate}{Resource replenishment rate}

\item{resource_capacity}{Resource carrying capacity}

\item{...}{Unused}
}
\value{
Vector containing resource spectrum at next timestep
}
\description{
If you set your resource dynamics to use this function then the time
evolution of the resource spectrum is described by a semi-chemostat equation
\deqn{\frac{\partial N_R(w,t)}{\partial t} = r_R(w) \Big[ c_R (w) - N_R(w,t) \Big] - \mu_R(w, t) N_R(w,t)}{dN_R(w,t)/d t  = r_R(w) ( c_R (w) - N_R(w,t) ) - \mu_R(w,t ) N_R(w,t)}
}
\details{
Here \eqn{r_R(w)} is the resource regeneration rate and \eqn{c_R(w)} is the
carrying capacity in the absence of predation. These parameters are changed
with \code{\link[=setResource]{setResource()}}. The mortality \eqn{\mu_R(w, t)} is
due to predation by consumers and is calculate with \code{\link[=getResourceMort]{getResourceMort()}}.

This function uses the analytic solution of the above equation to calculate
the resource abundance at time \code{t + dt} from all abundances and rates at time
\code{t}, keeping the mortality fixed during the timestep.

To set your model to use semichemostat dynamics for the resource you do

\if{html}{\out{<div class="sourceCode">}}\preformatted{params <- setResource(params, 
                      resource_dynamics = "resource_semichemostat", 
                      resource_level = 0.5)
}\if{html}{\out{</div>}}

where you should replace \code{params} with the name of the variable holding your
MizerParams object. You can of course choose any value between 0 and 1 for
the resource level.

The \code{\link[=balance_resource_semichemostat]{balance_resource_semichemostat()}} function is called by \code{\link[=setResource]{setResource()}}
to determine the values of the resource parameters that are needed to make
the replenishment rate at each size equal the consumption rate at that size,
as calculated by \code{\link[=getResourceMort]{getResourceMort()}}. It should be called with only one of
\code{resource_rate} or \code{resource_capacity} should and will return a named list
with the values for both.
}
\seealso{
Other resource dynamics: 
\code{\link{resource_constant}()},
\code{\link{resource_logistic}()}
}
\concept{resource dynamics}
