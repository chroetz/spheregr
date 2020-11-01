#' @param p 1x3
#' @param v kx3
#' @param phi nxk
#' @param y nx3
trigeo_energy <- function(p, v, phi, y) {
  z <- phi %*% v # nx3
  n <- nrow(z)
  norm_z <- sqrt(.rowSums(z^2, n, 3L, FALSE))
  co <- tcrossprod(y, p) * cos(norm_z) # nx1
  si <- .rowSums(y * z, n, 3L, FALSE) / norm_z * sin(norm_z) # nx1
  s <- clamp1(co + si)
  sum(acos(s)^2) / 2
}

trigeo_optim_fn <- function(par, phi, y, n_basis) {
  p <- convert_a2e_1(par[1:2])
  v <- matrix(NA, nrow=n_basis, ncol=3)
  for (j in 1:n_basis) v[j, ] <- get_v_rcpp(p, par[c(2*j+1, 2*j+2)])
  dim(p) <- c(1L, 3L)
  trigeo_energy(p, v, phi, y)
}

#' @param accuracy double in [0,1]. trade optimization accuracy for
#'   computational speed
trigeo_optim <- function(x, y, max_speed, n_basis, restarts = NULL, accuracy=0.3) {
  init_par <- get_initial_parameters(max_speed, restarts=restarts, n_basis=n_basis)
  res_lst <- list()
  phi <- t(eval_trig_base(x, n_basis+1))[,-1,drop=FALSE]
  for (i in seq_len(nrow(init_par))) {
    res_lst[[i]] <- stats::optim(
      init_par[i, ], trigeo_optim_fn,
      gr = NULL,
      phi = phi, y = y, n_basis = n_basis,
      method = "L-BFGS-B",
      lower = c(0, 0, rep(-max_speed, times=2*n_basis)),
      upper = c(pi, 2 * pi, rep(max_speed, times=2*n_basis)),
      control = list(factr = 10^(17 * (1 - accuracy)))
    )
  }
  values <- sapply(res_lst, function(x) x$value)
  idx <- which.min(values)
  res <- res_lst[[idx]]
  p <- convert_a2e_1(res$par[1:2])
  v <- matrix(NA, nrow=n_basis, ncol=3)
  for (j in 1:n_basis) v[j, ] <- get_v_rcpp(p, res$par[c(2*j+1, 2*j+2)])
  dim(p) <- c(1L, 3L)
  list(p = p, v = v, par=res$par)
}

#' @export
estimate_trigeo <- function(x, y, x_new, n_basis, periodize=FALSE, ...) {

  if (periodize) {
    n <- length(x)
    x <- c(x, rev(2-x))/2
    y <- rbind(y, y[n:1,])
    x_new <- x_new/2
  }

  res <- trigeo_optim(x, y, n_basis=n_basis, ...)
  phi <- t(eval_trig_base(x_new, n_basis+1))[,-1,drop=FALSE]
  estim <- Exp(res$p, phi %*% res$v)
  estim_a <- convert_e2a(estim)
  c(list(estim=estim, estim_a=estim_a), res)
}


#' @export
estimate_trigeo_loocv <- function(x, y, x_new, n_basis_max=5, ...) {
  n <- length(x)
  ns_basis <- 1:n_basis_max
  dists <- sapply(ns_basis, function(n_basis) {
    v <- sapply(seq_along(x), function(j) {
      res <- estimate_trigeo(x[-j], y[-j,], x[j], n_basis=n_basis, ...)
      dist(res$estim, y[j, ])
    })
    mean(v)
  })
  n_basis <- ns_basis[which.min(dists)]
  c(estimate_trigeo(x, y, x_new, n_basis=n_basis, ...), list(n_basis = n_basis))
}
