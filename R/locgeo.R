#' @param p 1x3
#' @param v 1x3
#' @param x nx1
#' @param y nx3
locgeo_energy <- function(p, v, x, y, w) {
  z <- x %*% v # nx3
  n <- nrow(z)
  norm_z <- sqrt(.rowSums(z^2, n, 3L, FALSE))
  co <- tcrossprod(y, p) * cos(norm_z) # nx1
  si <- .rowSums(y * z, n, 3L, FALSE) / norm_z * sin(norm_z) # nx1
  s <- clamp1(co + si)
  sum(acos(s)^2 * w) / 2
}

locgeo_optim_fn <- function(par, x, y, w) {
  p <- convert_a2e_1(par[1:2])
  v <- get_v_rcpp(p, par[3:4])
  dim(p) <- c(1L, 3L)
  dim(v) <- c(1L, 3L)
  locgeo_energy(p, v, x, y, w)
}

#' @param accuracy double in [0,1]. trade optimization accuracy for
#'   computational speed
locgeo_optim <- function(x, y, w, speed_bounds, restarts = 2, accuracy=0.3) {
  init_par <- get_initial_parameters(speed_bounds, restarts)
  res_lst <- list()
  for (i in seq_len(nrow(init_par))) {
    res_lst[[i]] <- stats::optim(
      init_par[i, ], locgeo_optim_fn,
      gr = NULL,
      x = x, y = y, w = w,
      method = "L-BFGS-B",
      lower = c(0, 0, -speed_bounds[2], -speed_bounds[2]),
      upper = c(pi, 2 * pi, speed_bounds[2], speed_bounds[2]),
      control = list(factr = 10^(17 * (1 - accuracy)))
    )
  }
  values <- sapply(res_lst, function(x) x$value)
  idx <- which.min(values)
  res <- res_lst[[idx]]
  p <- convert_a2e_1(res$par[1:2])
  v <- get_v_rcpp(p, res$par[3:4])
  dim(p) <- c(1L, 3L)
  dim(v) <- c(1L, 3L)
  list(p = p, v = v, par=res$par)
}

#' @export
estimate_locgeo <- function(x, y, x_new, kernel, h, speed_bounds) {
  estim <- matrix(NA_real_, ncol=3, nrow=length(x_new))
  for (i in seq_along(x_new)) {
    t <- x_new[i]
    w <- kernel((x-t)/h) / h
    w <- w / sum(w)
    if (!all(is.finite(w))) stop("Not all weights are finite.")
    res <- locgeo_optim(x, y, w, speed_bounds)
    estim[i,] <- Exp(res$p, t*res$v)
  }
  estim_a <- convert_e2a(estim)
  list(estim=estim, estim_a=estim_a)
}


#' @export
estimate_locgeo_loocv <- function(x, y, x_new, kernel, speed_bounds, n_h=7) {
  n <- length(x)
  hs <- (2/n)^seq(1, 0, len=n_h)
  dists <- sapply(hs, function(h) {
    v <- sapply(seq_along(x), function(j) {
      res <- estimate_locgeo(x[-j], y[-j,], x[j], kernel, h, speed_bounds)
      dist(res$estim, y[j, ])
    })
    mean(v)
  })
  h <- hs[which.min(dists)]
  c(estimate_locgeo(x, y, x_new, kernel, h, speed_bounds), list(h = h))
}
