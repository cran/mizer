## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE-----------------------------------------------------------
library(mizer)
params <- newTraitParams(no_sp = 2, min_w = 1e-3)
plotSpectra(params, resource = FALSE, power = 0)

