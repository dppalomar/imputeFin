#'
#' @title Fit Student's t VAR model to time series with missing values and/or outliers
#'
#' @description Estimate the parameters of a Student's t vector autoregressive
#'              model 
#'              \deqn{y_t = \phi_0 + \sum_{i=1}^p \Phi_i * y_{t-i} + \epsilon_t}
#'              to fit the given time series with missing values. 
#'              If the time series does not contain missing values, the 
#'              maximum likelihood (ML) estimation is done via the iterative
#'              EM algorithm until converge is achieved.
#'              With missing values, the stochastic EM algorithm is employed 
#'              for the estimation (currently the maximum number of iterations
#'              will be executed without attempting to check early converge).
#'
#' @param Y Time series object coercible to either a numeric matrix (e.g., \code{zoo} or \code{xts}) with missing values denoted by \code{NA}.
#' @param p Positive integer indicating the order of the VAR model.
#' @param omit_missing Logical value indicating whether to use the omit-variable method, i.e., 
#'                     excluding the variables with missing data from the analysis (default is \code{FALSE}).
#' @param parallel_max_cores Positive integer indicating the maximum numer of cores used in the parallel computing, 
#'                           only valid when \code{partition_groups} = \code{TRUE} (default is \eqn{1}).
#' @param verbose Logical value indicating whether to report in console the information of each iteration.
#' @param return_iterates Logical value indicating whether to return the parameter estimates at each iteration (default is \code{FALSE}).
#' @param initial List with the initial values of the parameters of the VAR model, which may contain some or all of the following elements:
#' \itemize{\item{\code{nu} (\eqn{\nu}) - a positive number as the degrees of freedom,}
#'          \item{\code{phi0} (\eqn{\phi_0}) - a numerical vector of length \code{ncol(Y)} as the interception of VAR model,}
#'          \item{\code{Phii} (\eqn{\Phi_i}) - a list of \code{p} matrices of dimension \code{ncol(Y)} as the autoregressive coefficient matrices,}
#'          \item{\code{scatter} (\eqn{\Sigma}) - a positive definite of dimension \code{ncol(Y)} as the scatter matrix.}}
#' @param L Positive integer with the number of Markov chains (default is \eqn{10}).
#' @param maxiter Positive integer with the number of maximum iterations (default is \eqn{100}).
#' @param ptol Non-negative number with the tolerance to determine the convergence of the (stochastic) EM method.
#' @param partition_groups Logical value indicating whether to partition \code{Y} into groups (default is \code{TRUE}).
#' @param K Positive integer indicating the values of the step sizes in the stochastic EM method.
#' 
#' @return A list with the following elements:
#' \item{\code{nu}}{The estimate for \eqn{\nu}.}
#' \item{\code{phi0}}{The estimate for \eqn{\phi_0}.}
#' \item{\code{Phii}}{The estimate for \eqn{\Phi_i}.}
#' \item{\code{scatter}}{The estimate for scatter matrix, i.e., \eqn{\Sigma}.}
#' 
#' \item{\code{converged}}{A logical value indicating whether the method has converged.}
#' \item{\code{iter_usage}}{A number indicating how many iteration has been used.}
#' \item{\code{elapsed_times}}{A numerical vector indicating how much is comsumed in each iteration.}
#' \item{\code{elapsed_time}}{A number indicating how much time is comsumed overall.}
#' \item{\code{elapsed_time_per_iter}}{A number indicating how much time is comsumed for each iteration in average.}
#' \item{\code{iterates_record}}{A list as the records of parameter estimates of each iteration, only returned when \code{return_iterates} = \code{TRUE}.}
#' 
#' @author Rui Zhou and Daniel P. Palomar
#' 
#' @references 
#' R. Zhou, J. Liu, S. Kumar, and D. P. Palomar, "Student’s t VAR Modeling with Missing Data via 
#' Stochastic EM and Gibbs Sampling," IEEE Trans. on Signal Processing, vol. 68, pp. 6198-6211, Oct. 2020.
#' 
#' @seealso \code{\link{fit_AR1_t}}
#' 
#' @examples 
#' library(imputeFin)
#' data(ts_VAR_t) 
#' fitted <- fit_VAR_t(Y = ts_VAR_t$Y, p = 2)
#' 
#' @importFrom mvtnorm dmvt
#' @importFrom magrittr %>% add
#' @importFrom utils head tail
#' @import parallel
#' @export
fit_VAR_t <- function(Y, p = 1, omit_missing = FALSE, parallel_max_cores = max(1, parallel::detectCores() - 1),
                      verbose = FALSE,
                      return_iterates = FALSE,
                      initial = NULL, L = 10, maxiter = 50, ptol = 1e-3, partition_groups = TRUE, K = round(maxiter/3)) {
  
  # error control
  if (!is.matrix(try(as.matrix(Y), silent = TRUE))) stop("\"Y\" must be coercible to a matrix.")
  if (p <= 0) stop("\"p\" must be a positive integer.")
  if (!is.logical(omit_missing)) stop("\"omit_missing\" must be a logical value.")
  if (parallel_max_cores <= 0) stop("\"parallel_max_cores\" must be a positive integer.")
  if (!is.logical(verbose)) stop("\"verbose\" must be a logical value.")
  if (!is.logical(return_iterates)) stop("\"return_iterates\" must be a logical value.")
  if (L <= 0) stop("\"L\" must be a positive integer.")
  if (maxiter < 1) stop("\"maxiter\" must be greater than 1.")
  if (ptol <= 0) stop("\"ptol\" must be greater than 0.")
  if (!is.logical(partition_groups)) stop("\"partition_groups\" must be a logical value.")
  if (K <= 0) stop("\"K\" must be a positive integer.")
  
  Y <- as.matrix(Y)
  T <- nrow(Y)
  N <- ncol(Y)
  # disassemble Y into Y_shreds, which is a list containing:
  # full_obs: a list of full observed segments 
  # part_obs: a list of partially observed segments
  if (partition_groups) {
    Y_shreds <- partitionMissingGroups(Y, p)
    if (omit_missing) {  # omit-variable method: ignore all partially observed segments
      Y_shreds$part_obs <- list()
      T <- sum(sapply(Y_shreds$full_obs, nrow) - p) + p  # NOTE: in this case, we need fake "T" since we have deleted some rows
    }
  } else {  # if choose not to partition into groups
    Y_shreds <- list(full_obs = list(), part_obs = list())
    if (anyNA(Y)) Y_shreds$part_obs <- list(Y) else Y_shreds$full_obs <- list(Y)
  }
  
  n_full_obs <- length(Y_shreds$full_obs)
  n_part_obs <- length(Y_shreds$part_obs)
  
  # decide number of parallel cores
  n_cores <- min(parallel_max_cores, parallel::detectCores() - 1, n_part_obs)
  # message(ls(environment() %>% parent.env(), all.names = TRUE))
  if (n_cores > 1) {
    if (verbose) message("creating a parallel socket cluster with ", n_cores, " cores...")
    cl <- parallel::makeCluster(n_cores)  # create parallel socket cluster
    clusterExport(cl = cl, varlist = assistant_fun_names, envir = environment() %>% parent.env())  # export assisting functions
    clusterEvalQ(cl = cl, expr = {library("magrittr")})  # dependencies
  } else {
    cl <- NULL
  }
  
  my_lapply <- function(DATA, FUN) {
    if (is.null(cl))
      lapply(DATA, FUN)
    else
      parallel::parLapplyLB(cl = cl, X = DATA, fun = FUN)
  }
  
  # browser()
  # initialize parameters
  nu      <- if (is.null(initial$nu)) 6 else initial$nu
  phi0    <- if (is.null(initial$phi0)) rep(0, N) else initial$phi0
  Phii    <- if (is.null(initial$Phii)) replicate(p, diag(N), FALSE) else initial$Phii
  scatter <- if (is.null(initial$scatter)) diag(N)  else initial$scatter
  
  # browser()
  Estep_part_obs_old <- EstepContainer(N, p)
  # chains_library <- lapply(Y_shreds$part_obs, function(part_obs) gibbsEstepMultiChains(part_obs, phi0, Phii, scatter, nu, L)$chain_states)  # init. markov chain states
  chains_library <- my_lapply(Y_shreds$part_obs, function(part_obs) gibbsEstepMultiChains(part_obs, phi0, Phii, scatter, nu, L)$chain_states)  # init. markov chain states
  
  # define tool functions for simplifying
  snapshot <- function() list(nu = nu, phi0 = phi0, Phii = Phii, scatter = scatter)
  s_tau    <- function() Estep$s_tau
  s_logtau <- function() Estep$s_logtau
  s_tauy   <- function(i) Estep$s_tauy[, i+1, drop = FALSE]
  S_tauyy  <- function(i, j) if (i >= j) Estep$s_tauyy[idx_mask(i, N), idx_mask(j, N)] else t(Estep$s_tauyy[idx_mask(j, N), idx_mask(i, N)])
  
  # let's loop
  if (return_iterates) iterates_record <- list(snapshot())
  elapsed_times <- c()
  
  for (iter in 1:maxiter) {
    if (verbose) {  # display the iteration info if required
      obj <- if (length(Y_shreds$part_obs) > 0) NA else my_lapply(Y_shreds$full_obs, function(full_obs_seg) loglikFullObs(full_obs_seg, phi0, Phii, scatter, nu)) %>% sum()
      message("iteration: ", iter, "\t objective:", obj) 
    }
    
    start_time <- proc.time()[3]  # record start time
    
    nu_old <- nu; phi0_old <- phi0; Phii_old <- Phii; scatter_old <- scatter  # record the current estimates
    
    # E-step ----------------------------------------
    # browser()
    # for the full observed segments
    if (n_full_obs > 0) {
      Estep_full_obs <- my_lapply(Y_shreds$full_obs, function(full_obs_seg) exactEstep(full_obs_seg, phi0, Phii, scatter, nu)) %>% do.call(add_num_list, .)
    } else {
      Estep_full_obs <- EstepContainer(N, p)
    }
    
    # for the partially observed segments
    if (n_part_obs > 0) {
      tmp <- my_lapply(1:n_part_obs, function(idx) gibbsEstepMultiChains(Y_shreds$part_obs[[idx]], phi0, Phii, scatter, nu, chains_library[[idx]]))
      chains_library <- lapply(tmp, function(x) x$chain_states)
      Estep_part_obs <- lapply(tmp, function(x) x$Estep) %>% do.call(add_num_list, .)
      
      # SAEM: combine the current Estep with the previous one
      convex_comb_para <- SAEMConvexCombPara(k = iter, K = K)
      if (convex_comb_para < 1) 
        Estep_part_obs <- add_num_list(mul_num_list(Estep_part_obs,     convex_comb_para),
                                        mul_num_list(Estep_part_obs_old, 1 - convex_comb_para))
      Estep_part_obs_old <- Estep_part_obs  # record the current Estep
    } else {
      Estep_part_obs <- EstepContainer(N, p)
    }
    Estep <- add_num_list(Estep_full_obs, Estep_part_obs)  
    
    
    # M-step ----------------------------------------
    # browser()
    # nu 
    Q_nu <- function(nu) ((nu/2)*log(nu/2) - log(gamma(nu/2)))*(T-p) + (nu/2) * (s_logtau() - s_tau())
    nu <- optimize(Q_nu, interval = c(2 + 1e-12, 100), maximum = TRUE)$maximum
    # if (!is.null(prior_info$nu)) nu <- prior_info$nu

    
    M <- assembleM(Estep = Estep)
    Psi <- M$M0 %*% solve(M$M1)
    phi0 <- Psi[, 1]
    Phii <- lapply(1:p, function(i) Psi[, idx_mask(i-1, N) + 1])
    
    S <- S_tauyy(0, 0) - 2 * M$M0 %*% t(Psi) + Psi %*% M$M1 %*% t(Psi)
    # browser()
    scatter <- S / (T-p)
    scatter <- (scatter + t(scatter)) / 2
    # Q(T, phi0, Phii, scatter, nu, Estep)
    
    if (return_iterates) iterates_record[[iter + 1]] <- snapshot()
    
    have_params_converged <-
      all(abs(fnu(nu) - fnu(nu_old))               <= .5 * ptol * (abs(fnu(nu_old)) + abs(fnu(nu)))) &&
      all(abs(phi0 - phi0_old)                     <= .5 * ptol * (abs(phi0) + abs(phi0_old))) &&
      all(abs(unlist(Phii) - unlist(Phii_old))     <= .5 * ptol * (abs(unlist(Phii)) + abs(unlist(Phii_old)))) &&
      all(abs(scatter - scatter_old)               <= .5 * ptol * (abs(scatter) + abs(scatter_old)))
    
    elapsed_times <- c(elapsed_times, proc.time()[3] - start_time)
    
    if (have_params_converged) break
  }
  
  if (!is.null(cl)) stopCluster(cl)
  
  # changed names of returned variables
  var_names <- colnames(Y)
  if (is.null(var_names)) var_names <- paste0("x", 1:N)
  for (i in 1:p) colnames(Phii[[i]]) <- rownames(Phii[[i]]) <- var_names
  colnames(scatter) <- rownames(scatter) <- var_names
  names(elapsed_times) <- NULL
  
  # return results -------------
  vars_to_be_returned <- list("nu"                    = nu, 
                              "phi0"                  = phi0, 
                              "Phii"                  = Phii, 
                              "scatter"               = scatter, 
                              "converged"             = (iter < maxiter),
                              "iter_usage"            = iter,
                              "elapsed_times"         = elapsed_times,
                              "elapsed_time"          = sum(elapsed_times),
                              "elapsed_time_per_iter" = mean(elapsed_times))
  if (return_iterates) {
    names(iterates_record) <- paste("iter", 0:(length(iterates_record)-1))
    vars_to_be_returned$iterates_record <- iterates_record
  }
  
  return(vars_to_be_returned)
  
}

fnu <- function(nu) nu/(nu-2)

SAEMConvexCombPara <- function(k, K) {
  if (k <= K) return(1)
  else return(1/(k - K))
}









###############################################################################
# assistant functions ---------
###############################################################################

assistant_fun_names <- 
  c("add_num_list", "assembleB", "assembleEstep", "assembleM" , "condGsnMoms", 
    "cutHeadNA", "EstepContainer", "exactEstep", "gibbsEstepMultiChains", 
    "gibbsEstepSingleChain", "gibbsTau", "gibbsY", "idx_mask", "lagInterception", 
    "mul_num_list", "partitionConsecutiveNonNA", 
    "partitionMissingGroups", ".Random.seed", "SAEMConvexCombPara", "sampleCondGsn")

# log-likelihood function when Y is full observed
loglikFullObs <- function(Y, phi0, Phii, scatter, nu) {
  dmvt(x = tail(Y, - length(Phii)) - lagInterception(Phii, Y, TRUE), delta = phi0, sigma = scatter, df = nu, log = TRUE) %>% sum()
}


# partition matrix Y according to missing patterns
# the sub-matrix with NAs should be identified and picked out
partitionMissingGroups <- function(Y, p) {
  T <- nrow(Y)
  obs_idx <- miss_idx <- rep(NA, T)
  
  # observation partition, when consecutive observations >= p + 1
  obs_mask <- !apply(Y, 1, anyNA)
  obs_idx[obs_mask] <- (1:T)[obs_mask]
  obs_partition <- partitionConsecutiveNonNA(x = obs_idx, consec = p + 1)
  
  # extract the raw missing mask
  miss_mask <- miss_mask_snapshort <- !obs_mask
  
  # amend the consecutive observations as missing (via using logical "TRUE" in the "miss_mask") when its length <= p
  tmp <- partitionConsecutiveNonNA(x = obs_idx, consec = 1)
  for (idx in tmp) if (length(idx) <= p) miss_mask[idx] <- TRUE
  
  # missing groups partition, include also p observations before and after
  miss_idx[miss_mask] <- (1:T)[miss_mask]
  miss_partition <- partitionConsecutiveNonNA(x = miss_idx, consec = 1)
  miss_partition <- lapply(miss_partition, function(x) {
    tmp <- c(min(x) + (-p:-1), x, max(x) + (1:p))
    tmp <- tmp[tmp <= T]
    tmp <- tmp[tmp >= 1]
    return(tmp)
    })
  
  # browser()
  # sanity check 
  # TODO: remove this part before submitting to CRAN 
  # tmp <- lapply(c(obs_partition, miss_partition), function(x) tail(x, -p))
  # tmp <- unlist(tmp)
  # if (length(tmp) != T-p || length(tmp) != length(unique(tmp))) stop("ERROR found in missing groups partition step!")
  
  return(list("full_obs" = lapply(obs_partition,  function(x) Y[x, ]), 
              "part_obs" = lapply(miss_partition, function(x) Y[x, ])))
}


# remove the heading NAs from a vector
# {examples} input:  c(NA, NA, 1, 2, 3)
#            output: c(1, 2, 3)
cutHeadNA <- function(x) {
  if (all(is.na(x))) {
    integer(0)
  } else if (!is.na(x[1])) {
    x
  } else {
    tail(x, 1-which(!is.na(x))[1] )
  }
}

# partition segments consisting of consecutive observations
# {example} input:  c(1, 2, 3, NA, 5, 6, 7, NA, NA, 10, 11, 12)
#           output: a list consisting of c(1, 2, 3), c(5, 6, 7), c(10, 11, 12)
partitionConsecutiveNonNA <- function(x, consec = 2) {
  x <- cutHeadNA(x)
  res <- list()
  count <- 1
  while (TRUE) {
    cutting <- which(is.na(x))[1] - 1  # the last position before the next missing 
    if (is.na(cutting)) cutting <- length(x)
    if (cutting >= consec) {  # found satisfying segment with length >= consec
      res[[count]] <- x[1:cutting]
      count <- count + 1
    }
    x <- tail(x, -cutting)  # drop the heading observations
    x <- cutHeadNA(x)  # drop the heading NAs
    if (length(x) == 0) break
  }
  return(res)
}

# the index mask function, particularly designed to start from 0
idx_mask <- function(i, N) (1:N) + i*N

# .loglikNoMiss <- function(X, phi0, Phii, delta, nu, scatter) {
#   # browser()
#   N <- ncol(X)
#   T <- nrow(X)
#   p <- length(Phii)
#   innovations <- tail(X, -p) - .crossprodPsiiX(Phii, X, TRUE) - matrix(data = phi0, nrow = T-p, ncol = N, byrow = TRUE)
#   sum(log(EMMIXuskew::dmst(dat = innovations, mu = rep(0, N), sigma = scatter, delta = delta, dof = nu)))
# }

# lagged interception from previous p observations 
# note the data matrix is ranked from old to recent
lagInterception <- function(Phii, Y, rm.last = FALSE) {
  T <- nrow(Y)
  p <- length(Phii)
  res <- do.call(cbind, lapply(1:p, function(idx) Y[(p+1-idx):(T+1-idx), , drop = FALSE])) %*% t(do.call(cbind, Phii))
  if (rm.last) head(res, -1) else res
}


# recursively add two numerical lists
# very useful when we collect the E-step statistics
# assuming the passed lists have the same structure (not check since it's an inner function)
add_num_list <- function(...) {
  data <- list(...)  # collect inputs
  res <- data[[1]][sapply(data[[1]], function(x) is.list(x) | is.numeric(x))]  # only allow numbers or lists
  if (length(data) == 1) return(res)
  for (i in 2:length(data)) {
    for (name in names(res)) {
      if (is.list(res[[name]]))
        res[[name]] <- add_num_list(res[[name]], data[[i]][[name]])
      if (is.numeric(res[[name]]))
        res[[name]] <- res[[name]] + data[[i]][[name]]
    }
  }
  return(res)
}

# recursively multiply a scale to a numerical list
mul_num_list <- function(lis, sca) {
  for (name in names(lis)) {
    if (is.list(lis[[name]]))
      lis[[name]] <- mul_num_list(lis[[name]], sca)
    if (is.numeric(lis[[name]]))
      lis[[name]] <- lis[[name]] * sca
  }
  return(lis)
}




# container for minimal sufficient statistics in E-step
# all elements are initialized as "0" when created
# note these values are already sum up over t
EstepContainer <- function(N, p) {
  list("s_logtau" = 0,
       "s_tau"    = 0,
       "s_tauy"   = matrix(0, N, p + 1),
       "s_tauyy"  = matrix(0, (p+1)*N, (p+1)*N))
}


# exact expectation computation when Y elements are all observed
exactEstep <- function(Y, phi0, Phii, scatter, nu) {
  if (anyNA(Y)) stop("Y can not contain any NA values!")
  if (nrow(Y) <= 1) stop("Y must have more than 2 rows!")
  
  T <- nrow(Y)
  N <- ncol(Y)
  p <- length(Phii)
  
  container <- EstepContainer(N, p)
  # browser()
  # as the tau follows the gamma distribution
  scatter_inv <- solve(scatter)
  tmp <- apply(t(tail(Y, -p)) - phi0 - t(lagInterception(Phii, Y, TRUE)), 2, function(x) as.numeric(x%*%scatter_inv%*%x))
  
  tau_vec_exp <- sapply(tmp, function(x) (nu+N)/(nu+x))
  logtau_vec_exp <- sapply(tmp, function(x) digamma((nu+N)/2) - log((nu+x)/2))
  return(assembleEstep(taus = tau_vec_exp, logtaus = logtau_vec_exp, Y = Y))
}


# to assemble B matrix from the current estimates,
# see equation (36) of the reference paper
assembleB <- function(phi0, Phii) {
  N <- length(phi0)
  p <- length(Phii)
  B <- matrix(0, N*p, N*p)
  for (i in 1:p) 
    B[idx_mask(0, N), idx_mask(i-1, N)] <- Phii[[i]]
  if (p > 1)
    diag(B[N + idx_mask(0, N*(p - 1)), idx_mask(0, N*(p - 1))]) <- 1
  
  return(B)
}


# for a segment, compute the gaussian moments of last (T-p) rows
# conditional on the first p rows, taus, and the current estimates
# see Lemma 2 of the reference paper
condGsnMoms <- function(phi0, Phii, scatter, Y_head_p, taus) {
  N <- length(phi0)
  p <- length(Phii)
  T_minus_p <- length(taus)
  if (nrow(Y_head_p) != p) stop("Invalid Phii or Y_head_p!")
  
  B <- assembleB(phi0, Phii)
  Bs <- list(B)
  for (i in 2:(max(T_minus_p, 2))) Bs[[i]] <- Bs[[i - 1]] %*% B
  powerB <- function(idx) if(idx == 0) diag(N*p) else Bs[[idx]]
  
  xp <- as.vector(t(Y_head_p[p:1, ]))
  
  mu <- sapply(0:(T_minus_p-1), function(x) powerB(x)[1:N, 1:N]%*%phi0) %>% apply(., 1, cumsum) %>% rbind() %>% t() %>%
    add(sapply(1:(T_minus_p), function(x) powerB(x)[1:N, ]%*%xp)) %>% as.vector()
  
  Sigma <- matrix(0, length(mu), length(mu))
  for (i in 1:T_minus_p) {
    for (j in 1:i) {
      tmp <- lapply(1:min(i, j), function(q) powerB(i-q)[1:N, 1:N]%*%scatter%*%t(powerB(j-q)[1:N, 1:N])/taus[q])
      Sigma[idx_mask(i-1, N), idx_mask(j-1, N)] <- Reduce('+', tmp)
    }
  }
  Sigma <- Sigma + t(Sigma)
  for (i in 1:T_minus_p) Sigma[idx_mask(i-1, N), idx_mask(i-1, N)] <- Sigma[idx_mask(i-1, N), idx_mask(i-1, N)]/2
  
  return(list("mu" = mu, "Sigma" = Sigma))
}


# when data is assumed to follow a jointly Gaussian distribution
# sample the missing data conditional on the observed part and the joint distribution parameter
# see Lemma 3 of the reference paper
#' @importFrom MASS mvrnorm
sampleCondGsn <- function(y, mu, Sigma) {
  if (!anyNA(y)) return(x)
  if (all(is.na(y))) stop("Invaild y!")
  
  missing_pattern <- is.na(y)
  mu_missing <- mu[missing_pattern] + Sigma[missing_pattern, !missing_pattern] %*% solve(Sigma[!missing_pattern, !missing_pattern]) %*% (y[!missing_pattern] - mu[!missing_pattern])
  Sigma_missing <- Sigma[missing_pattern, missing_pattern] - Sigma[missing_pattern, !missing_pattern] %*% solve(Sigma[!missing_pattern, !missing_pattern]) %*% Sigma[!missing_pattern, missing_pattern]
  y_missing <- MASS::mvrnorm(n = 1, mu = mu_missing, Sigma = Sigma_missing)
  y[missing_pattern] <- y_missing
  
  return(y)
}


# Gibbs sampling - taus
gibbsTau <- function(Y, phi0, Phii, scatter, nu, scatter_inv = NULL) {
  T <- nrow(Y)
  N <- ncol(Y)
  p <- length(Phii)
  if (is.null(scatter_inv)) scatter_inv <- solve(scatter)
  tmp <- apply(t(tail(Y, -p)) - phi0 - t(lagInterception(Phii, Y, rm.last = TRUE)), 2, function(x) as.numeric(x%*%scatter_inv%*%x))
  tau_vec <- sapply(tmp, function(x) rgamma(n = 1, shape = (nu+N)/2, rate = (x+nu)/2))
  return(tau_vec)
}

# Gibbs sampling - Y (missing part) in fashion of row by row
gibbsY <- function(tau_vec, Y, phi0, Phii, scatter, nu, chain_state) {
  
  if (!anyNA(Y)) return(Y)
  T <- nrow(Y)
  N <- ncol(Y)
  p <- length(Phii)
  
  # impute missing values in Y by single observation
  miss_idx <- (1:T)[apply(Y, 1, anyNA)]
  for (i in miss_idx) {  
    
    idx_first_p <- (-p:-1) + i
    idx_next_p  <- if (i == T) integer(0) else (i+1):min(i+p, T)
    # compute Gaussian moments
    gaus_mean_cov <- condGsnMoms(phi0, Phii, scatter, Y[idx_first_p, , drop = FALSE], tau_vec[c(i, idx_next_p) - p])
    # impute Y missing
    y <- sampleCondGsn(as.vector(t(rbind(Y[i, ], chain_state[idx_next_p, ]))), gaus_mean_cov$mu, gaus_mean_cov$Sigma)
    Y[i, ] <- head(y, N)
  }
  
  return(Y)
}


# assemble the minimal sufficient statistics from Gibbs samples
assembleEstep <- function(taus, logtaus = log(taus), Y) {
  N <- ncol(Y)
  T <- nrow(Y)
  p <- nrow(Y) - length(taus)  # p is infered
  idxes <- function(i, T, p) ((p+1):T) - i  # indexes for (T-p) rows, "i" is the offset
  # browser()
  container          <- EstepContainer(N, p)
  container$s_tau    <- sum(taus)
  container$s_logtau <- sum(logtaus)
  container$s_tauy   <- sapply(0:p, function(i) colSums(Y[idxes(i, T, p), , drop = FALSE] * taus))
  for (i in 0:p) {
    for (j in 0:p) {
      container$s_tauyy[idx_mask(i, N), idx_mask(j, N)] <- t(Y[idxes(i, T, p), , drop = FALSE] * taus) %*% Y[idxes(j, T, p), , drop = FALSE]
    }
  }
  
  # sanity check
  # TODO: remove this before pushing to CRAN
  # if(!isSymmetric.matrix(container$s_tauyy)) stop("The computed E(tau * y * y) is not symmetric!")
  
  return(container)
}


# gibbs sampling of taus and missing ys in single markov chain
# n_iter: number of iterations in the markov chain
# n_drop: number of iterations to be dropped at the begining
gibbsEstepSingleChain <- function(n_iter = 5e3, n_drop = 1e3, Y, phi0, Phii, scatter, nu, Y_iterator) {
  T <- nrow(Y)
  N <- ncol(Y)
  p <- length(Phii)
  
  # browser()
  container <- EstepContainer(N, p)
  
  # initialize
  if (missing(Y_iterator)) {  # Y_iterator is not given, meaning it is the initial of the markov chain
    tau_iterator <- rgamma(T-p, shape = nu/2, rate = nu/2)
    Y_iterator <- Y
    Y_iterator[is.na(Y_iterator)] <- mean(Y_iterator, na.rm = TRUE)  # simply initialize missing data in Y by mean value of observed ones
  }
  
  for (iter in 1:n_iter) {
    tau_iterator <- gibbsTau(Y_iterator, phi0, Phii, scatter, nu)
    Y_iterator   <- gibbsY(tau_iterator, Y, phi0, Phii, scatter, nu, Y_iterator)
    
    if (iter <= n_drop) next
    container <- add_num_list(container, mul_num_list(assembleEstep(taus = tau_iterator, Y = Y_iterator), 1 / (n_iter - n_drop)))
    
  }
  
  list("state"  = Y_iterator,
       "Estep"  = container)
}


# gibbs sampling of taus and missing ys in multiple markov chain
# chain_states: when is a number, it means number of markov chains, implying the initial status (initialization is going to be done in "gibbsEstepSingleChain()")
#               when is a list, it means the status of Y in markov chains

gibbsEstepMultiChains <- function(Y, phi0, Phii, scatter, nu, chain_states) { 
  # browser()
  if (is.numeric(chain_states))  # initialize chain states if the initial states are given in numerical
    tmp <- lapply(1:chain_states, function(x) gibbsEstepSingleChain(n_iter = 10, n_drop = 10, Y, phi0, Phii, scatter, nu))
  else                           # move multiple chains to next status
    tmp <- lapply(chain_states, function(chain_state) gibbsEstepSingleChain(n_iter = 1, n_drop = 0, Y, phi0, Phii, scatter, nu, chain_state))
  
  n_chain <- ifelse(is.numeric(chain_states), chain_states, length(chain_states))
  
  list("chain_states" = lapply(tmp, function(x) x$state),
       "Estep"        = lapply(tmp, function(x) x$Estep) %>% do.call(add_num_list, .) %>% mul_num_list(., 1/n_chain))
}






# to assemble M_0 and M_1 matrix from the minimal sufficient statistics,
# see equation (12) of the reference paper
# argument "Estep" is a directly returned result from "assembleEstep()"
assembleM <- function(Estep) {
  N <- nrow(Estep$s_tauy)
  p <- ncol(Estep$s_tauy) - 1
  M0 <- cbind(Estep$s_tauy[, 1], Estep$s_tauyy[idx_mask(0, N), -idx_mask(0, N)])
  M1 <- matrix(0, 1+N*p, 1+N*p)
  
  M1[1, 1]   <- Estep$s_tau
  M1[-1, 1]  <- as.vector(Estep$s_tauy[, -1])
  M1[1, -1]  <- as.vector(Estep$s_tauy[, -1])
  M1[-1, -1] <- Estep$s_tauyy[-idx_mask(0, N), -idx_mask(0, N)]
  
  return(list("M0" = M0, "M1" = M1))
}

