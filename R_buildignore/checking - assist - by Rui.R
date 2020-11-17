genParaVARp_asy <- function(N, p) {
  phi0    <<- rnorm(N)
  Phii    <<- replicate(p, cov(matrix(rnorm(N*N*2), 2*N, N)) / p / 2 + matrix(rnorm(N*N), N, N)/10, FALSE)
  scatter <<- cov(matrix(rnorm(N*N*2), 2*N, N))
  nu      <<- sample(4:7, 1)
}


genParaVARp <- function(N, p) {
  phi0    <<- rnorm(N)
  Phii    <<- replicate(p, cov(matrix(rnorm(N*N*2), 2*N, N)) / p / 2, FALSE)
  scatter <<- cov(matrix(rnorm(N*N*2), 2*N, N))
  nu      <<- sample(4:7, 1)
}


genDataVARp <- function(T, phi0, Phii, scatter, nu, pad = 10*T) {
  N <- length(phi0)
  p <- length(Phii)
  
  Y <- mvtnorm::rmvt(n = p, sigma = scatter, df = nu, delta = rep(0, N))
  residual <- mvtnorm::rmvt(n = T - p + pad, sigma = scatter, df = nu, delta = rep(0, N))
  # cat("log-likelihood of residuals is:", sum(dmvt(x = tail(residual, T-p), delta = rep(0, N), sigma = scatter, df = nu, log = TRUE)), "\n")
  # print(tail(residual, T-p))
  for (t in 1:nrow(residual)) {
    y_old <- tail(Y, p)
    y_t <- phi0 + as.vector(.lagInterception(Phii, y_old)) + residual[t, ]
    Y <- rbind(Y, y_t)
  }
  # residual <<- tail(residual, T)
  tmp <- tail(Y, T)
  colnames(tmp) <- paste0("x", 1:N)
  rownames(tmp) <- c()
  return(tmp)
}

corrupt <- function(Y, miss_pct = 0, miss_items = 0, out_pct = 0, out_items = 0, out_magntd = 10, p = 1) {
  T <- nrow(Y)
  N <- ncol(Y)
  n_corrupt_NA <- round(T * miss_pct)
  n_corrupt_OL <- round(T * out_pct)
  n_corrupt <- n_corrupt_NA + n_corrupt_OL
  mask_corrupt <- sample((p+1):T, n_corrupt)
  
  for (i in head(mask_corrupt, n_corrupt_NA))
    Y[i, sample(N, miss_items)] <- NA
  
  for (i in tail(mask_corrupt, n_corrupt_OL)) 
    Y[i, sample(N, out_items)] <- (rbinom(out_items, 1, 0.5) - 0.5) * 2 * out_magntd
  
  return(Y)
  
}


plotConvergence <- function(res_fit) {
  if (is.null(res_fit$iterates_record))
    stop("Fitting result does not contain iteration converge. Make sure to use \"return_iterates = TRUE\" when doing the fitting.")
  
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install package \"ggplot2\"", call. = FALSE)
  
  if (!requireNamespace("reshape2", quietly = TRUE))
    stop("Please install package \"reshape2\"", call. = FALSE)
  
  iteration <- nu <- nu_div_nu_2 <- value <- variable <- log_likelihood <- NULL  # hack to avoid problems with checks
  
  p_all <- list()
  data <- data.frame(iteration = 0:(length(res_fit$iterates_record)-1))
  if (!is.null(res_fit$iterates_record[[1]]$nu)) {
    data$nu <- sapply(res_fit$iterates_record, `[[`, "nu")
    p_nu <- ggplot2::ggplot(data, ggplot2::aes(x = iteration, y = nu)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::ggtitle("Convergence of nu")
    print(p_nu)
    p_all <- c(p_all, list(p_nu))
    
  }
  
  if (!is.null(res_fit$iterates_record[[1]]$phi0)) {
    mu_matrix <- sapply(res_fit$iterates_record, `[[`, "phi0")
    rownames(mu_matrix) <- paste0("phi0", 1:nrow(mu_matrix))
    data <- cbind(data, as.data.frame(t(mu_matrix)))
    # ggplot(data, aes(x = iteration, y = mu2)) +
    #   geom_line() + geom_point()
    p_mu <- ggplot2::ggplot(reshape2::melt(data, measure.vars = rownames(mu_matrix)), ggplot2::aes(x = iteration, y = value, col = variable)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::ggtitle("Convergence of phi0")
    print(p_mu)
    p_all <- c(p_all, list(p_mu))
  }
  
  if (!is.null(res_fit$iterates_record[[1]]$scatter)) {
    diag_scatter_matrix <- sapply(res_fit$iterates_record, function(x) diag(x$scatter))
    rownames(diag_scatter_matrix) <- paste0("scatter_", 1:nrow(diag_scatter_matrix))
    data <- cbind(data, as.data.frame(t(diag_scatter_matrix)))
    p_scatter <- ggplot2::ggplot(reshape2::melt(data, measure.vars = rownames(diag_scatter_matrix)), ggplot2::aes(x = iteration, y = value, col = variable)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::ggtitle("Convergence of scatter matrix")
    print(p_scatter)
    p_all <- c(p_all, list(p_scatter))
  }
  
  if (!is.null(res_fit$iterates_record[[1]]$Phii)) {
    for (i in 1:length(res_fit$iterates_record[[1]]$Phii)) {
      diag_scatter_matrix <- sapply(res_fit$iterates_record, function(x) diag(x$Phii[[i]]))
      rownames(diag_scatter_matrix) <- paste0("Phi", i, "_", 1:nrow(diag_scatter_matrix))
      data <- cbind(data, as.data.frame(t(diag_scatter_matrix)))
      p_scatter <- ggplot2::ggplot(reshape2::melt(data, measure.vars = rownames(diag_scatter_matrix)), ggplot2::aes(x = iteration, y = value, col = variable)) +
        ggplot2::geom_line() + ggplot2::geom_point() +
        ggplot2::ggtitle(paste0("Convergence of Phi-", i, " matrix"))
      print(p_scatter)
      p_all <- c(p_all, list(p_scatter))
    }
  }
  
  return(invisible(p_all))
}