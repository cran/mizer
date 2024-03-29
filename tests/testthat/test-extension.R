params <- NS_params
e <- globalenv()
e$test_dyn <- function(params, ...) {
    111
}
e$nt <- function(params, t, ...) {
    params@initial_n * t
}
e$resource_encounter <- function(params, n, n_pp, n_other, ...) {
    mizerEncounter(params, n = n, n_pp = n_other$resource, ...)
}
e$semichemostat <- function(params, n_other, rates, dt, component, ...) {
    c <- params@other_params[[component]]
    interaction <- params@species_params$interaction_resource
    mort <- as.vector(interaction  %*% rates$pred_rate)
    tmp <- c$rate * c$capacity / (c$rate + mort)
    return(tmp - (tmp - n_other[[component]]) * exp(-(c$rate + mort) * dt))
}

# setRateFunction works ----
test_that("setRateFunction works", {
    expect_error(setRateFunction(params, rate = "wrong", fun = "sum"),
               "The `rate` argument must be one of")
    xxx <- "xx"
    expect_error(setRateFunction(params, rate = "Mort", fun = "xxx"),
                 "There is no function ")
    p <- setRateFunction(params, "Mort", "mizerMort")
    expect_identical(params@rates_funcs, p@rates_funcs)
    p <- setRateFunction(params, rate = "EGrowth", fun = "test_dyn")
    expect_identical(p@rates_funcs[["EGrowth"]], "test_dyn")
    r <- mizerRates(p, n = initialN(p), 
                    n_pp = initialNResource(p), 
                    n_other = initialNOther(p),
                    effort = 0,
                    rates_fns = lapply(p@rates_funcs, get))
    expect_identical(r$e_growth, 111)
})

test_that("Time is passed correctly to rate functions", {
    params@rates_funcs$Encounter <- "nt"
    expect_identical(getEncounter(params, t = 2), nt(params, 2))
    params@rates_funcs$FeedingLevel <- "nt"
    expect_identical(getFeedingLevel(params, time_range = 2), nt(params, 2))
    
    gears <- unique(as.character(gear_params(params)$gear))
    effort <- array(0, dim = c(3, 4), 
                    dimnames = list(time = 2020:2022,
                                    gear = gears))
    sim <- project(params, effort = effort, dt = 1)
    expect_identical(getFeedingLevel(sim, time_range = 2021:2022)[1, , ],
                     nt(params, 2021))
    #TODO: extend this
})

# getRateFunction works ----
test_that("getRateFunction works", {
    expect_error(getRateFunction(params, "test"),
        "The `rate` argument must be one of")
    all <- getRateFunction(params)
    expect_type(all, "list")
    expect_named(all)
    expect_identical(all$Rates, "mizerRates")
    expect_identical(getRateFunction(params, rate = "Mort"), "mizerMort")
})

# other_params ----
test_that("We can set and get other params", {
    expect_length(other_params(params), 0)
    expect_error(other_params(params) <- 5,
                 "other_params should be a named list")
    expect_error(other_params(params) <- list(5),
                 "other_params should be a named list")
    other_params(params)$test <- 5
    expect_identical(other_params(params)$test, 5)
    other_params(params)$test <- NULL
    expect_length(other_params(params), 0)
    expect_null(other_params(params)$test)
})

# components ----
test_that("We can set, get and remove components", {
    expect_error(setComponent(params, "test", 1),
                 '"dynamics_fun" is missing')
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      encounter_fun = "test_dyn",
                      mort_fun = "test_dyn",
                      component_params = list(a = 2))
    expect_identical(p@other_dynamics, list(test = "test_dyn"))
    expect_identical(p@other_encounter, list(test = "test_dyn"))
    expect_identical(p@other_mort, list(test = "test_dyn"))
    expect_identical(p@other_params, list(test = list(a = 2)))
    comp <- getComponent(p, "test")
    expect_mapequal(comp, list(initial_value = 1,
                               encounter_fun = "test_dyn",
                               dynamics_fun = "test_dyn",
                               mort_fun = "test_dyn",
                               component_params = list(a = 2)))
    all <- getComponent(p)
    expect_type(all, "list")
    expect_identical(all$test, comp)
    p2 <- setComponent(p, "test2", 2, 
                      dynamics_fun = "test_dyn",
                      mort_fun = "test_dyn")
    all2 <- getComponent(p2)
    expect_length(all2, 2)
    expect_length(all2$test2, 5)
    expect_null(all2$test2$encounter_fun)
    p <- setComponent(p2, "test2", 1, "test_dyn", mort_fun = NULL)
    expect_null(getComponent(p, "test2")$mort_fun)
    expect_error(removeComponent(p2, "test3"),
                 "There is no component named test3")
    expect_null(getComponent(p2, "test3"))
    p1 <- removeComponent(p2, "test")
    d <- getComponent(p1, "test2")
    expect_length(p1@other_dynamics, 1)
    expect_length(p1@other_encounter, 0)
})

# initial values ----
test_that("We can set and get initial values for MizerParams", {
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn")
    expect_identical(initialNOther(p), list(test = 1))
    p <- setComponent(p, "test", list(a = 1, b = 2), 
                      dynamics_fun = "test_dyn")
    expect_identical(initialNOther(p), list(test = list(a = 1, b = 2)))
    initialNOther(p)$test <- 3
    expect_identical(initialNOther(p), list(test = 3))
    expect_error(initialNOther(p)$test2 <- 3,
                 "The following components do not exist: test2")
    p <- setComponent(p, "test2", 2, 
                      dynamics_fun = "test_dyn")
    expect_identical(initialNOther(p), list(test = 3, test2 = 2))
    expect_error(initialNOther(p) <- list(test = 4),
                 "Missing values for components test2")
    initialNOther(p)$test <- 4
    expect_identical(initialNOther(p)$test, 4)
    # test that we can get initial values from MizerSim object
    sim <- project(p, t_max = 0.2, t_save = 0.1)
    expect_identical(initialNOther(sim)$test, 4)
})

# encounter and mortality functions are called ----
test_that("encounter and mortality functions are called", {
    e <- getEncounter(params)
    m <- getMort(params)
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      encounter_fun = "test_dyn",
                      mort_fun = "test_dyn")
    expect_identical(getEncounter(p), e + 111)
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      mort_fun = "test_dyn")
    expect_identical(getMort(p), m + 111)
})

test_that("We can access simulation results", {
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      encounter_fun = "test_dyn",
                      mort_fun = "test_dyn")
    sim <- project(p, t_max = 0.2, t_save = 0.1)
    expect_identical(finalNOther(sim), list("test" = 111))
    expect_identical(NOther(sim)[2, ], list(111))
    expect_identical(NOther(sim)[1, ], list(1))
})

test_that("component can mimic resource", {
    params <- NS_params
    initialNResource(params) <- initialNResource(params) / 10
    component_params <- list(capacity = params@cc_pp,
                             rate = params@rr_pp)
    params2 <- params %>% 
        setComponent("resource",
                     initial_value = initialNResource(params),
                     dynamics_fun = "semichemostat",
                     component_params = component_params) %>% 
        setRateFunction("Encounter", "resource_encounter")
    sim <- project(params, t_max = 0.1, t_save = 0.1, effort = 0)
    sim2 <- project(params2, t_max = 0.1, t_save = 0.1, effort = 0)
    expect_identical(finalNResource(sim2), finalNOther(sim2)$resource)
    expect_identical(finalNResource(sim), finalNResource(sim2))
    expect_identical(finalN(sim), finalN(sim2))
})
