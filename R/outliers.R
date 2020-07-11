find_outliers_AR1_Gaussian <- function(y, fitted, outlier_prob_th = 0.001) {  # outlier_prob_th = 0.002 is too high and detects too many outliers
  index_obs <- which(!is.na(y))
  index_outliers <- NULL
  for (i in 2:length(index_obs)) {
    delta_i <- index_obs[i] - index_obs[i-1]
    mu_expected_i <- sum(fitted$phi1^( 0:(delta_i - 1) )) * fitted$phi0 + fitted$phi1^delta_i * y[index_obs[i-1]]
    delta_mu <- abs(y[index_obs[i]] - mu_expected_i)
    # compute the probability of that observation (or larger in magnitude)
    #prob_tail <- mvtnorm::pmvt(lower = -Inf, upper = -delta_mu, df = Inf, sigma = fitted$sigma2)
    #prob_tail <- stats::pt(-delta_mu/sqrt(fitted$sigma2), df = Inf)
    #prob_tail <- stats::pnorm(-delta_mu/sqrt(fitted$sigma2))
    prob_tail <- stats::pnorm(-delta_mu, sd = sqrt(fitted$sigma2))
    if (2*prob_tail < outlier_prob_th) {
      index_outliers <- c(index_outliers, index_obs[i])
      y[index_obs[i]] <- mu_expected_i  # this is a quick fix of the outlier so that it doesn't wrongly detect subsequent outliers
    }
  }
  return(index_outliers)
}



find_outliers_AR1_t <- function(y, fitted, outlier_prob_th = 0.001) {  # outlier_prob_th = 0.001 is too low and misses some outliers
  index_obs <- which(!is.na(y))
  index_outliers <- NULL
  for (i in 2:length(index_obs)) {
    delta_i <- index_obs[i] - index_obs[i-1]
    mu_expected_i <- sum(fitted$phi1^( 0:(delta_i - 1) )) * fitted$phi0 + fitted$phi1^delta_i * y[index_obs[i-1]]
    #upper <- max(y[index_obs[i]], 2*mu_expected_i - y[index_obs[i]])
    #lower <- min(y[index_obs[i]], 2*mu_expected_i - y[index_obs[i]])
    delta_mu <- abs(y[index_obs[i]] - mu_expected_i)
    # compute the probability of that observation (or larger in magnitude)
    #prob_tail <- mvtnorm::pmvt(lower = -Inf, upper = -delta_mu, df = max(round(fitted$nu), 1), sigma = fitted$sigma2)
    prob_tail <- stats::pt(-delta_mu/sqrt(fitted$sigma2), df = max(round(fitted$nu), 1))
    #prob_center <- mvtnorm::pmvt(lower = lower - mu_expected_i, upper = upper - mu_expected_i, df = max(round(fitted$nu), 1), sigma = fitted$sigma2)
    if (2*prob_tail < outlier_prob_th) {
      index_outliers <- c(index_outliers, index_obs[i])
      y[index_obs[i]] <- mu_expected_i  # this is a quick fix of the outlier so that it doesn't wrongly detect subsequent outliers
    }
  }
  return(index_outliers)
}