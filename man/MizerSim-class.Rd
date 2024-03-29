% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MizerSim-class.R
\docType{class}
\name{MizerSim-class}
\alias{MizerSim-class}
\title{A class to hold the results of a simulation}
\description{
A class that holds the results of projecting a \linkS4class{MizerParams}
object through time using \code{\link[=project]{project()}}.
}
\details{
A new \code{MizerSim} object can be created with the \code{\link[=MizerSim]{MizerSim()}}
constructor, but you will never have to do that because the object is
created automatically by \code{\link[=project]{project()}} when needed.

As a user you should never have to access the slots of a MizerSim object
directly. Instead there are a range of functions to extract the information.
\code{\link[=N]{N()}} and \code{\link[=NResource]{NResource()}} return arrays with the saved abundances of
the species and the resource population at size respectively. \code{\link[=getEffort]{getEffort()}}
returns the fishing effort of each gear through time.
\code{\link[=getTimes]{getTimes()}} returns the vector of times at which simulation results
were stored and \code{\link[=idxFinalT]{idxFinalT()}} returns the index with which to access
specifically the value at the final time in the arrays returned by the other
functions. \code{\link[=getParams]{getParams()}} returns the \code{MizerParams} object that was
passed to \code{project()}. There are also several
\link{summary_functions} and \link{plotting_functions}
available to explore the contents of a \code{MizerSim} object.

The arrays all have named dimensions. The names of the \code{time} dimension
denote the time in years. The names of the \code{w} dimension are weights in grams
rounded to three significant figures. The names of the \code{sp} dimension are the
same as the species name in the order specified in the species_params data
frame. The names of the \code{gear} dimension are the names of the gears, in the
same order as specified when setting up the \code{MizerParams} object.

Extensions of mizer can use the \code{n_other} slot to store the abundances of
other ecosystem components and these extensions should provide their own
functions for accessing that information.

The \code{MizerSim} class has changed since previous versions of mizer. To use
a \code{MizerSim} object created by a previous version, you need to upgrade it
with \code{\link[=upgradeSim]{upgradeSim()}}.
}
\section{Slots}{

\describe{
\item{\code{params}}{An object of type \linkS4class{MizerParams}.}

\item{\code{n}}{Three-dimensional array (time x species x size) that stores the
projected community number densities.}

\item{\code{n_pp}}{An array (time x size) that stores the projected resource number
densities.}

\item{\code{n_other}}{A list array (time x component) that stores the projected
values for other ecosystem components.}

\item{\code{effort}}{An array (time x gear) that stores the fishing effort by time
and gear.}
}}

