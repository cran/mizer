% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper_functions.R
\name{newTraitParams}
\alias{newTraitParams}
\title{Set up parameters for a trait-based multispecies model}
\usage{
newTraitParams(
  no_sp = 11,
  min_w_max = 10,
  max_w_max = 10^4,
  min_w = 10^(-3),
  max_w = max_w_max,
  eta = 10^(-0.6),
  min_w_mat = min_w_max * eta,
  no_w = round(log10(max_w_max/min_w) * 20 + 1),
  min_w_pp = 1e-10,
  w_pp_cutoff = min_w_mat,
  n = 2/3,
  p = n,
  lambda = 2.05,
  r_pp = 0.1,
  kappa = 0.005,
  alpha = 0.4,
  h = 40,
  beta = 100,
  sigma = 1.3,
  f0 = 0.6,
  fc = 0.25,
  ks = NA,
  gamma = NA,
  ext_mort_prop = 0,
  reproduction_level = 1/4,
  R_factor = deprecated(),
  gear_names = "knife_edge_gear",
  knife_edge_size = 1000,
  egg_size_scaling = FALSE,
  resource_scaling = FALSE,
  perfect_scaling = FALSE,
  min_w_inf = deprecated(),
  max_w_inf = deprecated()
)
}
\arguments{
\item{no_sp}{The number of species in the model.}

\item{min_w_max}{The maximum size of the smallest species in the
community. This will be rounded to lie on a grid point.}

\item{max_w_max}{The maximum size of the largest species in the community.
This will be rounded to lie on a grid point.}

\item{min_w}{The size of the the egg of the smallest species. This also
defines the start of the community size spectrum.}

\item{max_w}{The largest size in the model. By default this is set to the
largest maximum size \code{max_w_max}. Setting it to something larger
only makes sense if you plan to add larger species to the model later.}

\item{eta}{Ratio between maturity size and maximum size of a species.
Ignored if \code{min_w_mat} is supplied. Default is 10^(-0.6),
approximately 1/4.}

\item{min_w_mat}{The maturity size of the smallest species. Default value is
\code{eta * min_w_max}. This will be rounded to lie on a grid point.}

\item{no_w}{The number of size bins in the community spectrum. These bins
will be equally spaced on a logarithmic scale. Default value is such that
there are 20 bins for each factor of 10 in weight.}

\item{min_w_pp}{The smallest size of the resource spectrum. By default this
is set to the smallest value at which any of the consumers can feed.}

\item{w_pp_cutoff}{The largest size of the resource spectrum. Default value
is min_w_max unless \code{perfect_scaling = TRUE} when it is Inf.}

\item{n}{Scaling exponent of the maximum intake rate.}

\item{p}{Scaling exponent of the standard metabolic rate. By default this is
equal to the exponent \code{n}.}

\item{lambda}{Exponent of the abundance power law.}

\item{r_pp}{Growth rate parameter for the resource spectrum.}

\item{kappa}{Coefficient in abundance power law.}

\item{alpha}{The assimilation efficiency.}

\item{h}{Maximum food intake rate.}

\item{beta}{Preferred predator prey mass ratio.}

\item{sigma}{Width of prey size preference.}

\item{f0}{Expected average feeding level. Used to set \code{gamma}, the
coefficient in the search rate. Ignored if \code{gamma} is given
explicitly.}

\item{fc}{Critical feeding level. Used to determine \code{ks} if it is not given
explicitly.}

\item{ks}{Standard metabolism coefficient. If not provided, default will be
calculated from critical feeding level argument \code{fc}.}

\item{gamma}{Volumetric search rate. If not provided, default is determined
by \code{\link[=get_gamma_default]{get_gamma_default()}} using the value of \code{f0}.}

\item{ext_mort_prop}{The proportion of the total mortality that comes from
external mortality, i.e., from sources not explicitly modelled. A number in
the interval [0, 1).}

\item{reproduction_level}{A number between 0 an 1 that determines the
level of density dependence in reproduction, see \code{\link[=setBevertonHolt]{setBevertonHolt()}}.}

\item{R_factor}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use
\code{reproduction_level = 1 / R_factor} instead.}

\item{gear_names}{The names of the fishing gears for each species. A
character vector, the same length as the number of species.}

\item{knife_edge_size}{The minimum size at which the gear or gears select
fish. A single value for each gear or a vector with one value for each
gear.}

\item{egg_size_scaling}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
If TRUE, the egg size is a constant fraction of the
maximum size of each species. This fraction is \code{min_w / min_w_max}. If
FALSE, all species have the egg size \code{w_min}.}

\item{resource_scaling}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
If TRUE, the carrying capacity for larger resource
is reduced to compensate for the fact that fish eggs and larvae are
present in the same size range.}

\item{perfect_scaling}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
If TRUE then parameters are set so that the community
abundance, growth before reproduction and death are perfect power laws. In
particular all other scaling corrections are turned on.}

\item{min_w_inf}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} The argument has been
renamed to \code{min_w_max} to make it clearer that it refers to the maximum
size of a species not the von Bertalanffy asymptotic size parameter.}

\item{max_w_inf}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} The argument has been
renamed to \code{max_w_max}.}
}
\value{
An object of type \code{MizerParams}
}
\description{
This functions creates a \code{MizerParams} object describing a trait-based
model. This is a simplification of the general size-based model used in
\code{mizer} in which the species-specific parameters are the same for all
species, except for the maximum size, which is considered the most
important trait characterizing a species. Other parameters are related to the
maximum size. For example, the size at maturity is given by \code{w_max *
eta}, where \code{eta} is the same for all species. For the trait-based model
the number of species is not important. For applications of the trait-based
model see Andersen & Pedersen (2010). See the \code{mizer} website for more
details and examples of the trait-based model.
}
\details{
The function has many arguments, all of which have default values. Of
particular interest to the user are the number of species in the model and
the minimum and maximum sizes.

The characteristic weights of the smallest species are defined by
\code{min_w} (egg size), \code{min_w_mat} (maturity size) and
\code{min_w_max} (maximum size). The maximum sizes of
the \code{no_sp} species
are logarithmically evenly spaced, ranging from \code{min_w_max} to
\code{max_w_max}.
Similarly the maturity sizes of the species are logarithmically evenly
spaced, so that the ratio \code{eta} between maturity size and maximum
size is the same for all species. If \code{egg_size_scaling = TRUE} then also
the ratio between maximum size and egg size is the same for all species.
Otherwise all species have the same egg size.

In addition to setting up the parameters, this function also sets up an
initial condition that is close to steady state.

The search rate coefficient \code{gamma} is calculated using the expected
feeding level, \code{f0}.

The option of including fishing is given, but the steady state may loose its
natural stability if too much fishing is included. In such a case the user
may wish to include stabilising effects (like \code{reproduction_level}) to ensure
the steady state is stable. Fishing selectivity is modelled as a knife-edge
function with one parameter, \code{knife_edge_size}, which is the size at which
species are selected. Each species can either be fished by the same gear
(\code{knife_edge_size} has a length of 1) or by a different gear (the length of
\code{knife_edge_size} has the same length as the number of species and the order
of selectivity size is that of the maximum size).

The resulting \code{MizerParams} object can be projected forward using
\code{project()} like any other \code{MizerParams} object. When projecting
the model it may be necessary to reduce \code{dt} below 0.1 to avoid any
instabilities with the solver. You can check this by plotting the biomass or
abundance through time after the projection.
}
\examples{
params <- newTraitParams()
sim <- project(params, t_max = 5, effort = 0)
plotSpectra(sim)
}
\seealso{
Other functions for setting up models: 
\code{\link{newCommunityParams}()},
\code{\link{newMultispeciesParams}()},
\code{\link{newSingleSpeciesParams}()}
}
\concept{functions for setting up models}
