rm(list = ls())

source("R/fit_VAR_t.R")
source("R_buildignore/checking - assist - by Rui.R")

set.seed(3)
N <- 3
T <- 200
p <- 2
genParaVARp(N = N, p = p)
nu <- 5
Psi <- cbind(phi0, do.call(cbind, Phii))

Y <- genDataVARp(T = T, phi0 = phi0, Phii = Phii, scatter = scatter, nu = nu)
Y_crp <- corrupt(Y, miss_pct = 0.1, miss_items = 2, p = p)

Y <- xts::xts(Y_crp, seq(as.Date("2020-03-15"), length = T, by = "days"))


# create variable
ts_VAR_t <- list("Y"       = Y,
                 "phi0"    = phi0,
                 "Phii"    = Phii,
                 "scatter" = scatter,
                 "nu"      = nu)

# save the data to the package
save(ts_VAR_t, file = "data/ts_VAR_t.RData", version = 2, compress = "xz")
