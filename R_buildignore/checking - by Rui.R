# checking code: when Y is full observed ----------------------------------
rm(list = ls())

source("R/fit_t_VAR.R")
source("R_buildignore/checking - assist - by Rui.R")

# set.seed(3)
N <- 3
T <- 100
p <- 2
genParaVARp_asy(N = N, p = p)
nu <- 5
Psi <- cbind(phi0, do.call(cbind, Phii))

Y <- genDataVARp(T = T, phi0 = phi0, Phii = Phii, scatter = scatter, nu = nu)

fitted <- fit_VAR_t(Y = Y, p = p, return_iterates = TRUE)

plotConvergence(fitted)





# checking code: when Y is partially observed ----------------------------------
rm(list = ls())

source("R/fit_t_VAR.R")
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

fitted <- fit_VAR_t(Y = Y_crp, p = p, max_iter = 50, return_iterates = TRUE)
fitted$elapsed_time_per_iter


plotConvergence(fitted)


# checking code: simple omit-variable method ----------------------------------
rm(list = ls())

source("R/fit_t_VAR.R")
source("R_buildignore/checking - assist - by Rui.R")

# set.seed(3)
N <- 3
T <- 200
p <- 2
genParaVARp(N = N, p = p)
nu <- 5
Psi <- cbind(phi0, do.call(cbind, Phii))

Y <- genDataVARp(T = T, phi0 = phi0, Phii = Phii, scatter = scatter, nu = nu)
Y_crp <- corrupt(Y, miss_pct = 0.1, miss_items = 2, p = p)

fitted <- fit_VAR_t(Y = Y_crp, p = p, omit_missing = TRUE, return_iterates = TRUE)

plotConvergence(fitted)




# checking code: when Y is partially observed (parallel mode) ----------------------------------
rm(list = ls())

source("R/fit_t_VAR.R")
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

fitted <- fit_VAR_t(Y = Y_crp, p = p, max_iter = 100, parallel_max_cores = 5, return_iterates = TRUE)
fitted$elapsed_time_per_iter

plotConvergence(fitted)






# .partitionMissingGroups -----------
Y <- matrix(1:20, 20, 1)
Y[4:5] <- Y[8] <- Y[13:16] <- Y[18] <- NA
Y
.partitionMissingGroups(Y, 1)
.partitionMissingGroups(Y, 2)
.partitionMissingGroups(Y, 3)


# .lagInterception ------------------