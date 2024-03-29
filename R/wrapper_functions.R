# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#' Set up parameters for a community-type model
#' 
#' This functions creates a \code{\linkS4class{MizerParams}} object describing a
#' community-type model. 
#' The function has many arguments, all of which have default values.
#' 
#' A community model has several features that distinguish it from a
#' multi-species model:
#' 
#' * Species identities of individuals are ignored. All are aggregated into a
#'   single community.
#' * The resource spectrum only extends to the start of the community spectrum. 
#' * Reproductive rate is constant, independent of the energy invested in 
#'   reproduction, which is set to 0. 
#' * Standard metabolism is turned off (the parameter `ks` is set to 0).
#'   Consequently, the growth rate is now determined solely by the assimilated
#'   food 
#' 
#' Fishing selectivity is modelled as a knife-edge function with one parameter, 
#' `knife_edge_size`, which determines the size at which species are 
#' selected.
#' 
#' The resulting `MizerParams` object can be projected forward using 
#' \code{project()} like any other `MizerParams` object. When projecting 
#' the community model it may be necessary to keep a small time step size
#' `dt` of around 0.1 to avoid any instabilities with the solver. You can
#' check for these numerical instabilities by plotting the biomass or abundance
#' through time after the projection.
#' 
#' @param max_w The maximum size of the community. The `w_max` of the 
#'   species used to represent the community is set to this value.
#' @param min_w The minimum size of the community.
#' @param z0 The background mortality of the community.
#' @param alpha The assimilation efficiency of the community.
#' @param f0 The average feeding level of individuals who feed on a power-law 
#'   spectrum. This value is used to calculate the search rate parameter 
#'   `gamma`.
#' @param h The coefficient of the maximum food intake rate.
#' @param n The allometric growth exponent. Used as allometric exponent for
#'   the maximum intake rate of the community as well as the intrinsic growth
#'   rate of the resource.
#' @param beta The preferred predator prey mass ratio.
#' @param sigma The width of the prey preference.
#' @param gamma Volumetric search rate. Estimated using `h`, `f0` and 
#'   `kappa` if not supplied.
#' @param reproduction The constant reproduction in the smallest size class of
#'   the community spectrum. By default this is set so that the community
#'   spectrum is continuous with the resource spectrum.
#' @param knife_edge_size The size at the edge of the knife-edge-selectivity
#'   function.
#' @inheritParams newMultispeciesParams
#' @export
#' @return An object of type \code{\linkS4class{MizerParams}}
#' @references K. H. Andersen,J. E. Beyer and P. Lundberg, 2009, Trophic and 
#'   individual efficiencies of size-structured communities, Proceedings of the 
#'   Royal Society, 276, 109-114
#' @family functions for setting up models
#' @examples
#' params <- newCommunityParams()
#' sim <- project(params, t_max = 10)
#' plotBiomass(sim)
#' plotSpectra(sim, power = 2)
#' 
#' # More satiation. More mortality
#' params <- newCommunityParams(f0 = 0.8, z0 = 0.4)
#' sim <- project(params, t_max = 10)
#' plotBiomass(sim)
#' plotSpectra(sim, power = 2)
newCommunityParams <- function(max_w = 1e6,
                               min_w = 1e-3,
                               no_w = 100,
                               min_w_pp = 1e-10,
                               z0 = 0.1,
                               alpha = 0.2,
                               f0 = 0.7,
                               h = 10,
                               gamma = NA,
                               beta = 100,
                               sigma = 2.0,
                               n = 2 / 3,
                               kappa = 1000,
                               lambda = 2.05,
                               r_pp = 10,
                               knife_edge_size = 1000,
                               reproduction) {
    w_max <- max_w
    w_pp_cutoff <- min_w
    ks <- 0 # Turn off standard metabolism
    p <- n # But not used as ks = 0
    
    # Make the species data.frame
    species_params <- data.frame(
        species = "Community",
        w_max = w_max,
        f0 = f0,
        h = h, # max food intake
        gamma = gamma, # vol. search rate,
        ks = ks, # standard metabolism coefficient,
        beta = beta,
        sigma = sigma,
        z0 = z0, # background mortality
        alpha = alpha,
        erepro = 1, # not used
        interaction_resource = 1,
        sel_func = "knife_edge",
        knife_edge_size = knife_edge_size,
        stringsAsFactors = FALSE
    )
    params <- 
        newMultispeciesParams(species_params, no_w = no_w, min_w_pp = min_w_pp,
                              p = p, n = n, lambda = lambda, min_w = min_w,
                              resource_capacity = kappa,
                              resource_rate = r_pp,
                              w_pp_cutoff = w_pp_cutoff,
                              info_level = 0)
    
    initial_n <- array(kappa * params@w ^ (-lambda), 
                       dim = c(1, length(params@w)))
    initialN(params) <- initial_n

    params@rates_funcs$RDD <- "constantRDD"
    if (missing(reproduction)) {
        reproduction <- get_required_reproduction(params)
    }
    params@species_params$constant_reproduction <- reproduction
    params@psi[] <- 0 # Need to force to be 0. Can try setting w_mat but 
                          # due to slope still not 0
    comment(params@psi) <- "No investment into reproduction in community model."
    # Set w_mat to NA for clarity - it is not actually being used
    params@species_params$w_mat[] <- NA
    return(params)
}

#' Set up parameters for a trait-based multispecies model
#' 
#' This functions creates a `MizerParams` object describing a trait-based
#' model. This is a simplification of the general size-based model used in
#' `mizer` in which the species-specific parameters are the same for all
#' species, except for the maximum size, which is considered the most
#' important trait characterizing a species. Other parameters are related to the
#' maximum size. For example, the size at maturity is given by \code{w_max *
#' eta}, where `eta` is the same for all species. For the trait-based model
#' the number of species is not important. For applications of the trait-based
#' model see Andersen & Pedersen (2010). See the `mizer` website for more
#' details and examples of the trait-based model.
#'
#' The function has many arguments, all of which have default values. Of
#' particular interest to the user are the number of species in the model and
#' the minimum and maximum sizes.
#'
#' The characteristic weights of the smallest species are defined by
#' `min_w` (egg size), `min_w_mat` (maturity size) and 
#' `min_w_max` (maximum size). The maximum sizes of 
#' the `no_sp` species
#' are logarithmically evenly spaced, ranging from `min_w_max` to
#' `max_w_max`. 
#' Similarly the maturity sizes of the species are logarithmically evenly
#' spaced, so that the ratio `eta` between maturity size and maximum
#' size is the same for all species. If \code{egg_size_scaling = TRUE} then also
#' the ratio between maximum size and egg size is the same for all species.
#' Otherwise all species have the same egg size.
#'
#' In addition to setting up the parameters, this function also sets up an
#' initial condition that is close to steady state.
#'
#' The search rate coefficient `gamma` is calculated using the expected
#' feeding level, `f0`.
#'
#' The option of including fishing is given, but the steady state may loose its
#' natural stability if too much fishing is included. In such a case the user
#' may wish to include stabilising effects (like `reproduction_level`) to ensure
#' the steady state is stable. Fishing selectivity is modelled as a knife-edge
#' function with one parameter, `knife_edge_size`, which is the size at which
#' species are selected. Each species can either be fished by the same gear
#' (`knife_edge_size` has a length of 1) or by a different gear (the length of
#' `knife_edge_size` has the same length as the number of species and the order
#' of selectivity size is that of the maximum size).
#'
#' The resulting `MizerParams` object can be projected forward using
#' \code{project()} like any other `MizerParams` object. When projecting
#' the model it may be necessary to reduce `dt` below 0.1 to avoid any
#' instabilities with the solver. You can check this by plotting the biomass or
#' abundance through time after the projection.
#'
#' @param no_sp The number of species in the model.
#' @param min_w_max The maximum size of the smallest species in the
#'   community. This will be rounded to lie on a grid point.
#' @param max_w_max The maximum size of the largest species in the community.
#'   This will be rounded to lie on a grid point.
#' @param min_w The size of the the egg of the smallest species. This also
#'   defines the start of the community size spectrum.
#' @param max_w The largest size in the model. By default this is set to the
#'   largest maximum size `max_w_max`. Setting it to something larger
#'   only makes sense if you plan to add larger species to the model later.
#' @param eta Ratio between maturity size and maximum size of a species.
#'   Ignored if `min_w_mat` is supplied. Default is 10^(-0.6),
#'   approximately 1/4.
#' @param min_w_mat The maturity size of the smallest species. Default value is
#'   \code{eta * min_w_max}. This will be rounded to lie on a grid point.
#' @param no_w The number of size bins in the community spectrum. These bins
#'   will be equally spaced on a logarithmic scale. Default value is such that
#'   there are 20 bins for each factor of 10 in weight.
#' @param min_w_pp The smallest size of the resource spectrum. By default this
#'   is set to the smallest value at which any of the consumers can feed.
#' @param w_pp_cutoff The largest size of the resource spectrum. Default value
#'   is min_w_max unless \code{perfect_scaling = TRUE} when it is Inf.
#' @param n Scaling exponent of the maximum intake rate.
#' @param p Scaling exponent of the standard metabolic rate. By default this is
#'   equal to the exponent `n`.
#' @param lambda Exponent of the abundance power law.
#' @param r_pp Growth rate parameter for the resource spectrum.
#' @param kappa Coefficient in abundance power law.
#' @param alpha The assimilation efficiency.
#' @param ks Standard metabolism coefficient. If not provided, default will be
#'   calculated from critical feeding level argument `fc`.
#' @param fc Critical feeding level. Used to determine `ks` if it is not given
#'   explicitly.
#' @param h Maximum food intake rate.
#' @param beta Preferred predator prey mass ratio.
#' @param sigma Width of prey size preference.
#' @param f0 Expected average feeding level. Used to set `gamma`, the
#'   coefficient in the search rate. Ignored if `gamma` is given
#'   explicitly.
#' @param gamma Volumetric search rate. If not provided, default is determined
#'   by [get_gamma_default()] using the value of `f0`.
#' @param ext_mort_prop The proportion of the total mortality that comes from
#'   external mortality, i.e., from sources not explicitly modelled. A number in
#'   the interval [0, 1).
#' @param reproduction_level A number between 0 an 1 that determines the
#'   level of density dependence in reproduction, see [setBevertonHolt()].
#' @param R_factor `r lifecycle::badge("deprecated")` Use
#'   `reproduction_level = 1 / R_factor` instead.
#' @param gear_names The names of the fishing gears for each species. A
#'   character vector, the same length as the number of species.
#' @param knife_edge_size The minimum size at which the gear or gears select
#'   fish. A single value for each gear or a vector with one value for each
#'   gear.
#' @param egg_size_scaling `r lifecycle::badge("experimental")`
#'   If TRUE, the egg size is a constant fraction of the
#'   maximum size of each species. This fraction is \code{min_w / min_w_max}. If
#'   FALSE, all species have the egg size `w_min`.
#' @param resource_scaling `r lifecycle::badge("experimental")`
#'   If TRUE, the carrying capacity for larger resource
#'   is reduced to compensate for the fact that fish eggs and larvae are
#'   present in the same size range.
#' @param perfect_scaling `r lifecycle::badge("experimental")`
#'   If TRUE then parameters are set so that the community
#'   abundance, growth before reproduction and death are perfect power laws. In
#'   particular all other scaling corrections are turned on. 
#' @param min_w_inf `r lifecycle::badge("deprecated")` The argument has been
#'   renamed to `min_w_max` to make it clearer that it refers to the maximum
#'   size of a species not the von Bertalanffy asymptotic size parameter.
#' @param max_w_inf `r lifecycle::badge("deprecated")` The argument has been
#'   renamed to `max_w_max`.
#' @export
#' @importFrom lifecycle deprecated
#' @return An object of type `MizerParams`
#' @family functions for setting up models
#' @examples
#' params <- newTraitParams()
#' sim <- project(params, t_max = 5, effort = 0)
#' plotSpectra(sim)
newTraitParams <- function(no_sp = 11,
                           min_w_max = 10,
                           max_w_max = 10 ^ 4,
                           min_w = 10 ^ (-3),
                           max_w = max_w_max,
                           eta = 10^(-0.6),
                           min_w_mat = min_w_max * eta,
                           no_w = round(log10(max_w_max / min_w) * 20 + 1),
                           min_w_pp = 1e-10,
                           w_pp_cutoff = min_w_mat,
                           n = 2 / 3,
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
                           reproduction_level = 1 / 4,
                           R_factor = deprecated(),
                           gear_names = "knife_edge_gear",
                           knife_edge_size = 1000,
                           egg_size_scaling = FALSE,
                           resource_scaling = FALSE,
                           perfect_scaling = FALSE,
                           min_w_inf = deprecated(),
                           max_w_inf = deprecated()) {
    ## Deprecated arguments ----
    if (lifecycle::is_present(min_w_inf)) {
        lifecycle::deprecate_warn(
            when = "2.4.0.0", 
            what = "newTraitParams(min_w_inf)",
            with = "newTraitParams(min_w_max)"
        )
        min_w_max <- min_w_inf
    }
    if (lifecycle::is_present(max_w_inf)) {
        lifecycle::deprecate_warn(
            when = "2.4.0.0", 
            what = "newTraitParams(max_w_inf)",
            with = "newTraitParams(max_w_max)"
        )
        max_w_max <- max_w_inf
    }
    
    ## Check validity of parameters ----
    assert_that(is.flag(egg_size_scaling),
                is.flag(resource_scaling),
                is.flag(perfect_scaling))
    if (ext_mort_prop >= 1 || ext_mort_prop < 0) {
        stop("ext_mort_prop must be a number between 0 and 1",
             " because it should be the proportion of the total mortality",
             " coming from sources other than predation.")
    }
    if (min_w <= 0) {
        stop("The smallest egg size min_w must be greater than zero.")
    }
    no_w <- round(no_w)
    if (no_w < log10(max_w_max / min_w)*5) {
        no_w <- round(log10(max_w_max / min_w) * 5 + 1)
        message(paste("Increased no_w to", no_w, "so that there are 5 bins ",
                      "for an interval from w and 10w."))
    }
    if (no_w > 10000) {
        message("Running a simulation with ", no_w, 
                " size bins is going to be very slow.")
    }
    if (min_w_max >= max_w_max) {
        stop("The maximum size of the smallest species min_w_max must be ",
             "smaller than the maximum size of the largest species max_w_max")
    }
    if (min_w >= min_w_mat) {
        stop("The egg size of the smallest species min_w must be smaller than ",
             "its maturity size min_w_mat")
    }
    if (min_w_mat >= min_w_max) {
        stop("The maturity size of the smallest species min_w_mat must be ",
             "smaller than its maximum size min_w_max")
    }
    no_sp <- as.integer(no_sp)
    if (no_sp < 2) {
        stop("The number of species must be at least 2.")
    }
    if (!all(c(n, r_pp, lambda, kappa, alpha, h, beta, sigma, f0) > 0)) {
        stop("The parameters n, lambda, r_pp, kappa, alpha, h, beta, sigma ",
             "and f0, if supplied, need to be positive.")
    }
    if (!is.na(fc) && (fc < 0 || fc > f0)) {
        stop("The critical feeding level must lie between 0 and f0")
    }
    if (!is.na(gamma)) {  # If gamma is supplied, f0 is ignored
        f0 <- NA
    }
    # Check gears
    if (length(knife_edge_size) > no_sp) {
        stop("knife_edge_size needs to be no longer than the number of species in the model")
    }
    if ((length(knife_edge_size) > 1) && (length(knife_edge_size) != no_sp)) {
        warning("Length of knife_edge_size is less than number of species so gear information is being recycled. Is this what you want?")
    }
    if ((length(gear_names) != 1) && (length(gear_names) != no_sp)) {
        stop("Length of gear_names argument must equal the number of species.")
    }
    if (lifecycle::is_present(R_factor)) {
        lifecycle::deprecate_soft("2.3.0", "newTraitParams(R_factor)",
                                  "newTraitParams(reproduction_level)",
                                  "Set `reproduction_level = 1 / R_factor`.")
        reproduction_level <- 1 / R_factor
    }
    
    if (perfect_scaling) {
        egg_size_scaling <- TRUE
        resource_scaling <- TRUE
        w_pp_cutoff <- Inf
        p <- n
    }
    
    ## Set grid points and characteristic sizes ----
    # in such a way that the sizes all line up with the grid and the species are
    # all equally spaced.
    
    # Divide the range from min_w to max_w into (no_w - 1) logarithmic bins of
    # log size dx so that the last bin starts at max_w
    min_x <- log10(min_w)
    max_x <- log10(max_w)
    dx <- (max_x - min_x) / (no_w - 1) 
    x <- seq(min_x, by = dx, length.out = no_w)
    w <- 10 ^ x
    
    # Find index of nearest grid point to min_w_max that is an integer multiple
    # of the species spacing away from max_w
    min_x_inf <- log10(min_w_max)
    max_x_inf <- log10(max_w_max)
    # bins_per_sp is the number of bins separating species
    bins_per_sp <- round((max_x_inf - min_x_inf) / (dx * (no_sp - 1)))
    min_i_inf <- no_w - (no_sp - 1) * bins_per_sp
    # Maximum sizes for all species
    w_max <- w[seq(min_i_inf, by = bins_per_sp, length.out = no_sp)]
    
    # Find index of nearest grid point to min_w_mat
    min_x_mat <- log10(min_w_mat)
    min_i_mat <- round((min_x_mat - min_x) / dx) + 1
    # Maturity sizes for all species
    w_mat <- w[seq(min_i_mat, by = bins_per_sp, length.out = no_sp)]
    
    if (egg_size_scaling) {
        # Determine egg weights w_min for all species
        w_min_idx <- seq(1, by = bins_per_sp, length.out = no_sp)
        w_min <- w[w_min_idx]
    } else {
        w_min <- rep(min_w, no_sp)
        w_min_idx <- rep(1, no_sp)
    }

    ## Build Params Object ----
    erepro <- 0.1  # Will be changed later to achieve coexistence
    species_params <- data.frame(
        species = as.factor(1:no_sp),
        w_min = w_min,
        w_max = w_max,
        w_mat = w_mat,
        w_min_idx = w_min_idx,
        h = h,
        gamma = gamma,
        ks = ks,
        f0 = f0,
        fc = fc,
        beta = beta,
        sigma = sigma,
        z0 = 0,
        alpha = alpha,
        erepro = erepro,
        stringsAsFactors = FALSE
    )
    gear_params <- data.frame(
        gear = gear_names,
        species = species_params$species,
        sel_func = "knife_edge",
        knife_edge_size = knife_edge_size,
        catchability = 1,
        stringsAsFactors = FALSE
    )
    params <-
        newMultispeciesParams(
            species_params,
            gear_params = gear_params,
            min_w = min_w,
            no_w = no_w,
            max_w = max_w,
            w_pp_cutoff = max_w,
            lambda = lambda,
            kappa = kappa,
            n = n,
            p = p,
            min_w_pp = min_w_pp,
            resource_rate = r_pp,
            info_level = 0
        )
    
    w <- params@w
    dw <- params@dw
    w_full <- params@w_full
    ks <- params@species_params$ks[[1]]
    f0 <- get_f0_default(params)[[1]]
    
    ## Construct steady state solution ----
    
    # Get constants for steady-state solution
    # Predation mortality rate coefficient
    mu0 <- get_power_law_mort(params)
    # Add backgound mortality rate
    mu0 <- mu0 / (1 - ext_mort_prop)
    hbar <- alpha * h * f0 - ks
    if (hbar < 0) {
        stop("The feeding level is not sufficient to maintain the fish.")
    }
    pow <- mu0 / hbar / (1 - n)
    if (pow < 1) {
        message("The ratio of death rate to growth rate is too small, leading to
                an accumulation of fish at their largest size.")
    }
    
    initial_n <- params@psi  # get array with correct dimensions and names
    initial_n[, ] <- 0
    mumu <- mu0 * w^(n - 1)  # Death rate
    i_inf <- min_i_inf  # index of maximum size
    i_min <- 1  # index of natural egg size
    for (i in 1:no_sp) {
        gg <- hbar * w^n * (1 - params@psi[i, ])  # Growth rate
        idx <- w_min_idx[i]:(i_inf - 2)
        # Steady state solution of the upwind-difference scheme used in project
        n_exact <- c(1, cumprod(gg[idx] / ((gg + mumu * dw)[idx + 1])))
        # Use the first species for normalisation
        if (i == 1) {
            dist_sp <- bins_per_sp * dx
            mult <- kappa / 
                sum(n_exact * (w ^ (lambda - 1) * dw)[1:(min_i_inf - 1)]) *
                (10 ^ (dist_sp * (1 - lambda) / 2) - 10^(-dist_sp * (1 - lambda) / 2)) / 
                (1 - lambda)
        }
        if (!egg_size_scaling) {
            n_exact <- n_exact / n_exact[i_min]
        }
        idxs <- w_min_idx[i]:(i_inf - 1) 
        initial_n[i, idxs] <- n_exact * mult * 
            (w_max[1] / w_max[i]) ^ lambda
        i_inf <- i_inf + bins_per_sp
        i_min <- i_min + bins_per_sp
    }
    params@initial_n <- initial_n
    
    # Calculate the community spectrum
    sc <- colSums(initial_n)
    params@sc <- sc
    
    ##  Setup resource ----
    if (resource_scaling) {
        resource_vec <- (kappa * w ^ (-lambda)) - sc
        # Cut off resource at w_pp_cutoff
        resource_vec[w >= w_pp_cutoff] <- 0
        if (any(resource_vec < 0)) {
            if (!perfect_scaling) {
                # Do not allow negative resource abundance
                message("Note: Negative resource abundance values overwritten with zeros")
                resource_vec[resource_vec < 0] <- 0
            } else {
                message("Note: Negative resource abundances")
            }
        }
        params@cc_pp[sum(params@w_full <= w[1]):length(params@cc_pp)] <-
            resource_vec
    }
    if (!perfect_scaling) {
        params@cc_pp[w_full >= w_pp_cutoff] <- 0
    }
    params@initial_n_pp <- params@cc_pp
    
    # The cc_pp factor needs to be higher than the desired steady state in
    # order to compensate for predation mortality
    m2_background <- getResourceMort(params)
    params@cc_pp <- (params@rr_pp + m2_background ) * 
        params@initial_n_pp / params@rr_pp
    comment(params@cc_pp) <- paste("Carrying capacity set higher than the",
                                   "desired steady state in order to",
                                   "compensate for predation mortality.")
    
    ## Setup external death ----
    m2 <- getPredMort(params)
    flag <- FALSE
    for (i in 1:no_sp) {
        # The steplike psi was only needed when we wanted to use the analytic
        # expression for the steady-state solution
        # params@psi[i,] <- (w / w_max[i]) ^ (1 - n)
        # params@psi[i, w < (w_mat[i] - 1e-10)] <- 0
        # params@psi[i, w > (w_max[i] - 1e-10)] <- 1
        params@mu_b[i,] <- mu0 * w ^ (n - 1) - m2[i, ]
        if (!perfect_scaling && any(params@mu_b[i,] < 0)) {
            params@mu_b[i, params@mu_b[i,] < 0] <- 0
        }
        comment(params@mu_b) <- paste("Background mortality is set so that in",
                                      "combination with predation mortality",
                                      "it gives a power-law mortality.")
    }
    
    
    ## Set reproduction to meet boundary condition ----
    params@species_params$erepro <- params@species_params$erepro *
        get_required_reproduction(params) / getRDI(params) 
    params <- setBevertonHolt(params, reproduction_level = reproduction_level)

    return(params)
}


# Helper function to calculate the coefficient of the death rate created by
# a power-law spectrum of predators, assuming they have the same predation 
# parameters as the first species.
get_power_law_mort <- function(params) {
    params@interaction[] <- 0
    params@interaction[1, 1] <- 1
    params@initial_n[1, ] <- params@resource_params$kappa * 
        params@w^(-params@resource_params$lambda)
    return(getPredMort(params)[1, 1] / 
               params@w[[1]] ^ (1 + params@species_params[[1, "q"]] - 
                                    params@resource_params$lambda))
}

#' Determine reproduction rate needed for initial egg abundance
#'
#' @param params A MizerParams object
#' @return A vector of reproduction rates for all species
#' @concept helper
get_required_reproduction <- function(params) {
    assert_that(is(params, "MizerParams"))
    
    no_sp <- nrow(params@species_params)
    mumu <- getMort(params)
    gg <- getEGrowth(params)
    reproduction <- params@species_params$erepro # vector of correct length
    for (i in (1:no_sp)) {
        gg0 <- gg[i, params@w_min_idx[i]]
        mumu0 <- mumu[i, params@w_min_idx[i]]
        DW <- params@dw[params@w_min_idx[i]]
        reproduction[i] <- params@initial_n[i, params@w_min_idx[i]] *
            (gg0 + DW * mumu0)
    }
    return(reproduction)
}

