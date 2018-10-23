#
# Univariate estimation functions
#
EstimateMoments_uv_iid <- function(x) {
  tmp <- Estimate_AR1(x, phi1=0)
  return(list(mu=tmp$phi0, var=tmp$var, nu=tmp$nu))
}

EstimateMomentsResidual_uv_randomwalk <- function(x) {
  tmp <- Estimate_AR1(x, phi1=1)
  return(list(mu=tmp$phi0, var=tmp$var, nu=tmp$nu))
}

Estimate_AR1 <- function(x, phi0=NULL, phi1=NULL, var=NULL, nu=NULL) {
  x <- as.matrix(x)
  estimate_phi0 <- ifelse(is.null(phi0), TRUE, FALSE)
  estimate_phi1 <- ifelse(is.null(phi1), TRUE, FALSE)
  estimate_nu <- ifelse(is.null(nu), TRUE, FALSE)
  estimate_var <- ifelse(is.null(var), TRUE, FALSE)

  #
  # code here...
  #
    
  return( list(phi0=phi0, phi1=phi1, var=var, nu=nu) )
}

EstimateMomentsResidual_ARp <- function() {
  
}


#
# Multivariate estimation functions
#

EstimateMoments_mv_iid <- function() {
  
}

EstimateMomentsResidual_mv_randomwalk <- function() {
  
}

Estimate_VAR1 <- function() {
  
}

Estimate_VARp <- function() {
  
}