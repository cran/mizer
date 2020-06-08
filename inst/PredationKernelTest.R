rm(list=ls())

library(devtools)
install_github("sizespectrum/mizer", ref = "dev")
library(mizer)

mariaParams <- read.csv(file = "inst/Tasm_Params.csv") #16 species
inter <- read.csv(file = "inst/Tasm_Inter.csv") 
inter <- as.matrix(inter)

## setup basic parameters 
no_size_groups = 200
tmax = 100
dt = 0.2
kappa = 0.02 # 1.5 #1.5 ## 0.017g/m3 
lambda = 2.1 
w_pp_cutoff = 10 #5
r_pp = 0.4
min_w_pp = 1e-10

#setup standard params
params <- MizerParams(mariaParams, interaction = inter, no_w = no_size_groups, kappa = kappa, lambda = lambda, w_pp_cutoff = w_pp_cutoff, r_pp = r_pp, min_w_pp = min_w_pp)

model <- project(params, t_max = tmax, effort= 0, dt = dt)
plot(model)
getBiomass(model)[tmax,]

## make new predation kernel as per help example  
beta <- params@species_params$beta
sigma <- params@species_params$sigma
w <- params@w
w_full <- params@w_full
pk = array(beta, dim = c(length(beta),length(w), length(w_full)))
pk <- exp(-0.5 * sweep(log(sweep(sweep(pk, 3, w_full, "*") ^ -1, 2, w, "*")),
                       1, sigma, "/") ^ 2)
pk <- sweep(pk, c(2, 3), combn(w_full, 1, function(x, w) x < w, w = w), "*")

### try with new predation kernel 

params1 <- change_pred_kernel(params, pk)
model1 <- project(params1, t_max = tmax, effort= 0, dt = dt)
plot(model1)
getBiomass(model1)[tmax,]

## Just look at the biomass ratio at the last step 
round((getBiomass(model)[tmax,]/getBiomass(model1)[tmax,]),3)

no_size_groups = 1000
params <- MizerParams(mariaParams, interaction = inter, no_w = no_size_groups,
                      kappa = kappa, lambda = lambda, w_pp_cutoff = w_pp_cutoff,
                      r_pp = r_pp, min_w_pp = min_w_pp)

model_1000 <- project(params, t_max = tmax, effort= 0, dt = dt)

w <- params@w
w_full <- params@w_full
pk = array(beta, dim = c(length(beta),length(w), length(w_full)))
pk <- exp(-0.5 * sweep(log(sweep(sweep(pk, 3, w_full, "*") ^ -1, 2, w, "*")),
                       1, sigma, "/") ^ 2)
pk <- sweep(pk, c(2, 3), combn(w_full, 1, function(x, w) x < w, w = w), "*")

### try with new predation kernel 

params1 <- change_pred_kernel(params, pk)
model1_1000 <- project(params1, t_max = tmax, effort= 0, dt = dt)

round((getBiomass(model1_1000)[tmax,]/getBiomass(model1)[tmax,]),3)
