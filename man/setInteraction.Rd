% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setInteraction.R
\name{setInteraction}
\alias{setInteraction}
\alias{interaction_matrix}
\alias{interaction_matrix<-}
\title{Set species interaction matrix}
\usage{
setInteraction(params, interaction = NULL)

interaction_matrix(params)

interaction_matrix(params) <- value
}
\arguments{
\item{params}{MizerParams object}

\item{interaction}{Optional interaction matrix of the species (predator
species x prey species). By default all entries are 1. See "Setting
interaction matrix" section below.}

\item{value}{An interaction matrix}
}
\value{
\code{setInteraction}: A MizerParams object with updated interaction
matrix

\code{interaction_matrix()}: The interaction matrix (predator species x
prey species)
}
\description{
Set species interaction matrix
}
\section{Setting interaction matrix}{


You do not need to specify an interaction matrix. If you do not, then the
predator-prey interactions are purely determined by the size of predator
and prey and totally independent of the species of predator and prey.

The interaction matrix \eqn{\theta_{ij}} modifies the interaction of each
pair of species in the model. This can be used for example to allow for
different spatial overlap among the species.
The values in the interaction matrix are used to scale the encountered food
and predation mortality (see on the website \href{https://sizespectrum.org/mizer/articles/model_description.html#sec:pref}{the section on predator-prey encounter rate}
and on \href{https://sizespectrum.org/mizer/articles/model_description.html#mortality}{predation mortality}).
The first index refers to the predator species and the second to the prey
species.

The interaction matrix is used when calculating the food encounter rate in
\code{\link[=getEncounter]{getEncounter()}} and the predation mortality rate in \code{\link[=getPredMort]{getPredMort()}}. Its
entries are dimensionless numbers. If all the values in the interaction
matrix are equal then predator-prey interactions are determined entirely by
size-preference.

This function checks that the supplied interaction matrix is valid and then
stores it in the \code{interaction} slot of the \code{params} object.

The order of the columns and rows of the \code{interaction} argument should be
the same as the order in the species params data frame in the \code{params}
object. If you supply a named array then the function will check the order
and warn if it is different. One way of creating your own interaction
matrix is to enter the data using a spreadsheet program and saving it as a
.csv file. The data can then be read into R using the command \code{read.csv()}.

The interaction of the species with the resource are set via a column
\code{interaction_resource} in the \code{species_params} data frame. By default this
column is set to all 1s.
}

\examples{
params <- newTraitParams(no_sp = 3)
inter <- getInteraction(params)
inter[1, 2:3] <- 0
params <- setInteraction(params, interaction = inter)
getInteraction(params)
}
\seealso{
Other functions for setting parameters: 
\code{\link{gear_params}()},
\code{\link{setExtEncounter}()},
\code{\link{setExtMort}()},
\code{\link{setFishing}()},
\code{\link{setInitialValues}()},
\code{\link{setMaxIntakeRate}()},
\code{\link{setMetabolicRate}()},
\code{\link{setParams}()},
\code{\link{setPredKernel}()},
\code{\link{setReproduction}()},
\code{\link{setSearchVolume}()},
\code{\link{species_params}()}
}
\concept{functions for setting parameters}
