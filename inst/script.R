version <- "2.0"

pc <- set_community_model(max_w = 1e7, min_w = 0.0001, z0 = 0.2, alpha = 0.1,
                          h = 12, beta = 50, sigma = 1.9, q = 0.85, n = 0.7,
                          kappa = 2000, f0 = 0.6, r_pp = 1,
                          knife_edge_size = 500)
simc <- project(pc, effort = 0.5, t_max = 10, dt = 0.05, t_save = 5)
plotSpectra(simc)
saveRDS(simc, paste0("simc.", version, ".rds"))

pcs <- set_community_model()
simcs <- project(pcs)
plotSpectra(simcs)
saveRDS(simcs, paste0("simcs.", version, ".rds"))

pt <- set_trait_model(no_sp = 3, min_w_inf = 15, max_w_inf = 1e6,
                      no_w = 50, min_w = 0.002, min_w_pp = 1e-9,
                      w_pp_cutoff = 0.1, k0 = 40, n = 0.7, p = 0.8, 
                      q = 0.8, eta = 0.2, r_pp = 10, kappa = 0.001,
                      alpha = 0.2, ks = 5, z0pre = 0.5, h = 40,
                      beta = 60, sigma = 1.5, f0 = 0.65,
                      knife_edge_size = 500)
simt <- project(pt, effort = 0.5, t_max = 10, dt = 0.05, t_save = 5)
plotSpectra(simt)
saveRDS(simt, paste0("simt.", version, ".rds"))

pts <- set_trait_model()
simts <- project(pts)
plotSpectra(simts)
saveRDS(simts, paste0("simts.", version, ".rds"))

data("NS_species_params")
data(inter)
pm <- MizerParams(NS_species_params, inter, max_w = 1e7, min_w = 0.0001,
                  no_w = 50, min_w_pp = 1e-9, n = 0.7, p = 0.75, q = 0.85,
                  r_pp = 5, kappa = 1e10, w_pp_cutoff = 1, f0 = 0.8,
                  z0pre = 0.5, z0exp = 0.2)
simm <- project(pm)
plotSpectra(simm)
saveRDS(simm, paste0("simm.", version, ".rds"))

pms <- MizerParams(NS_species_params)
simms <- project(pms, t_max = 10, t_Save = 5)
plotSpectra(simms)
saveRDS(simms, paste0("simms.", version, ".rds"))

sim <- upgradeSim(simm.2.0)
validObject(sim)
sim1 <- upgradeSim(simm.1.0)
validObject(sim1)
diffobj::diffObj(sim1@params, sim@params)
diffobj::diffObj(sim1, sim)

p1 <- plotSpectra(sim1)
p2 <- plotSpectra(sim)
vdiffr::widget_slide(p1, p2)
vdiffr::widget_toggle(p1, p2)

params <- sim@params
params1 <- sim1@params
simn <- project(params, effort = 0.5, t_max = 10, dt = 0.05, t_save = 5)
simn1 <- project(params1, effort = 0.5, t_max = 10, dt = 0.05, t_save = 5)
p1 <- plotSpectra(simn)
p2 <- plotSpectra(simn1)
vdiffr::widget_slide(p1, p2)
