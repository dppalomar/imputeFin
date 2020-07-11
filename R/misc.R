
# extracts the diagonal on top of the main diagonal
diag1 <- function(X) {
  m <- min(dim(X))
  X[1 + dim(X)[1L] + 0L:(m - 2L) * (dim(X)[1L] + 1)]  # main diag: x[1 + 0L:(m - 1L) * (dim(x)[1L] + 1)]
}


is_inner_NA <- function(y) {
  if (NCOL(y) > 1) stop("Function is_inner_NA() cannot deal with multiple columns.")
  is_na_y <- is.na(y)
  index_obs  <- which(!is_na_y)
  index_obs_min <- min(index_obs)
  index_obs_max <- max(index_obs)
  is_na_y[1:index_obs_min] <- FALSE
  is_na_y[index_obs_max:length(y)] <- FALSE
  return(is_na_y)
}


any_inner_NA <- function(y) {
  idx_obs <- which(!is.na(y))
  if (length(idx_obs) == 0)  # all NAs!
    return(FALSE)
  else
    anyNA(y[min(idx_obs):max(idx_obs)])
}

