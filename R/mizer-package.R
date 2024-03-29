#' mizer: Multi-species size-based modelling in R
#'
#' The mizer package implements multi-species size-based modelling in R. It has 
#' been designed for modelling marine ecosystems.
#'
#' Using \pkg{mizer} is relatively simple.  There are three main stages: 
#' 
#' 1. *Setting the model parameters*. This is done by creating an object of
#' class \linkS4class{MizerParams}. This includes model parameters such as the
#' life history parameters of each species, and the range of the size spectrum.
#' There are several setup functions that help to create a MizerParams objects
#' for particular types of models:
#'     * [newSingleSpeciesParams()]
#'     * [newCommunityParams()]
#'     * [newTraitParams()]
#'     * [newMultispeciesParams()]
#' 
#' 2. *Running a simulation*. This is done by calling the
#' [project()] function with the model parameters. This produces an
#' object of \linkS4class{MizerSim} that contains the results of the simulation.
#'
#' 3. *Exploring results*. After a simulation has been run, the results can be
#' explored using a range of [plotting_functions], [summary_functions] and
#' [indicator_functions].
#'
#' See the [mizer website](https://sizespectrum.org/mizer/) for full details of
#' the principles behind mizer and how the package can be used to perform
#' size-based modelling.
#'
#' @import ggplot2 methods assertthat dplyr
#' @importFrom plotly ggplotly plotlyOutput renderPlotly
#' @importFrom stats fft mvfft lm pnorm runif complete.cases
#' @importFrom grDevices col2rgb
#' @importFrom utils modifyList packageVersion globalVariables
#' @importFrom rlang signal cnd_muffle
#' @importFrom lifecycle deprecated
"_PACKAGE"

#' @importFrom reshape2 melt
#' @export
reshape2::melt

globalVariables(c("expect_equal", "observed", "model", "is_observed"))
