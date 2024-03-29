#' Set species interaction matrix
#'
#' @section Setting interaction matrix:
#'
#'   You do not need to specify an interaction matrix. If you do not, then the
#'   predator-prey interactions are purely determined by the size of predator
#'   and prey and totally independent of the species of predator and prey.
#'
#'   The interaction matrix \eqn{\theta_{ij}} modifies the interaction of each
#'   pair of species in the model. This can be used for example to allow for
#'   different spatial overlap among the species.
#'   The values in the interaction matrix are used to scale the encountered food
#'   and predation mortality (see on the website [the section on predator-prey
#'   encounter
#'   rate](https://sizespectrum.org/mizer/articles/model_description.html#sec:pref)
#'    and on [predation
#'   mortality](https://sizespectrum.org/mizer/articles/model_description.html#mortality)).
#'    The first index refers to the predator species and the second to the prey
#'   species.
#'
#'   The interaction matrix is used when calculating the food encounter rate in
#'   [getEncounter()] and the predation mortality rate in [getPredMort()]. Its
#'   entries are dimensionless numbers. If all the values in the interaction
#'   matrix are equal then predator-prey interactions are determined entirely by
#'   size-preference.
#'
#'   This function checks that the supplied interaction matrix is valid and then
#'   stores it in the `interaction` slot of the `params` object.
#'
#'   The order of the columns and rows of the `interaction` argument should be
#'   the same as the order in the species params data frame in the `params`
#'   object. If you supply a named array then the function will check the order
#'   and warn if it is different. One way of creating your own interaction
#'   matrix is to enter the data using a spreadsheet program and saving it as a
#'   .csv file. The data can then be read into R using the command `read.csv()`.
#'
#'   The interaction of the species with the resource are set via a column
#'   `interaction_resource` in the `species_params` data frame. By default this
#'   column is set to all 1s.
#'
#' @param params MizerParams object
#' @param interaction Optional interaction matrix of the species (predator
#'   species x prey species). By default all entries are 1. See "Setting
#'   interaction matrix" section below.
#'
#' @return `setInteraction`: A MizerParams object with updated interaction
#'   matrix
#' @export
#' @family functions for setting parameters
#' @examples
#' params <- newTraitParams(no_sp = 3)
#' inter <- getInteraction(params)
#' inter[1, 2:3] <- 0
#' params <- setInteraction(params, interaction = inter)
#' getInteraction(params)
setInteraction <- function(params,
                           interaction = NULL) {
    assert_that(is(params, "MizerParams"))
    if (is.null(interaction)) {
        interaction <- params@interaction
    }
    if (!is.matrix(interaction)) {
        interaction <- as.matrix(interaction)
    }
    if (!is.numeric(interaction)) {
        stop("The entries of the interaction matrix should be numeric.")
    }
    # Check dims of interaction argument
    if (!identical(dim(params@interaction), dim(interaction))) {
        stop("interaction matrix is not of the right dimensions. Must be ",
             "number of species x number of species.")
    }
    # Check that all values of interaction matrix are positive
    if (!all(interaction >= 0)) {
        stop("All entries in the interaction matrix must be non-negative.")
    }
    # In case user has supplied names to interaction matrix, check them.
    if (!is.null(dimnames(interaction))) {
        if (!is.null(names(dimnames(interaction)))) {
            if (!identical(names(dimnames(interaction)),
                           names(dimnames(params@interaction)))) {
                message("Note: Your interaction matrix has dimensions called: `",
                        toString(names(dimnames(interaction))),
                        "`. I expected 'predator, prey'. ", 
                        "I will now ignore your names.")
            }
        }
        names(dimnames(interaction)) <- names(dimnames(params@interaction))
        # If user did not supply rownames, then save to assume that they have
        # put the rows in the same order as the columns, so copy over 
        # the colnames
        if (is.null(rownames(interaction)) || 
            all(rownames(interaction) == 
                as.character(seq_len(nrow(interaction))))) {
            rownames(interaction) <- colnames(interaction)
        }
        if (!identical(dimnames(params@interaction)[[1]],
                       dimnames(interaction)[[1]]) ||
            !identical(make.names(dimnames(params@interaction)[[2]]),
                       make.names(dimnames(interaction)[[2]]))) {
            message("Note: Dimnames of interaction matrix do not match the ",
                    "order of species names in the species data.frame. I am ",
                    "now ignoring your dimnames so your interaction matrix ",
                    "may be in the wrong order.")
        }
    }
    params@interaction[] <- interaction
    
    # Check the interaction_resource column in species_params
    message <- "Note: No interaction_resource column in species data frame so assuming all species feed on resource."
    species_params <- set_species_param_default(params@species_params,
                                                "interaction_resource", 1,
                                                message = message)
    # Check that all values of interaction vector are positive
    if (!all(species_params$interaction_resource >= 0)) {
        stop("Values in the resource interaction vector must be non-negative.")
    }
    params@species_params$interaction_resource <- species_params$interaction_resource
    
    params@time_modified <- lubridate::now()
    return(params)
}

#' Deprecated function to get interaction matrix
#' 
#' You should now use [interaction_matrix()] instead.
#' 
#' @param params A MizerParams object
#' @export
#' @keywords internal
getInteraction <- function(params) {
    lifecycle::deprecate_warn("2.4.0", "getInteraction()", 
                              "interaction_matrix()")
    interaction_matrix(params)
}


#' @rdname setInteraction
#' @return `interaction_matrix()`: The interaction matrix (predator species x
#'   prey species)
#' @export
interaction_matrix <- function(params) {
    params@interaction
}

#' @rdname setInteraction
#' @param value An interaction matrix
#' @export
`interaction_matrix<-` <- function(params, value) {
    setInteraction(params, interaction = value)
}
