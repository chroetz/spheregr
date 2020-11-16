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

locgeo_optim <- function(x, y, w, max_speed, grid_size, accuracy) {
  init_par <- get_initial_parameters(grid_size, max_speed)
  res_lst <- list()
  for (i in seq_len(nrow(init_par))) {
    res_lst[[i]] <- stats::optim(
      init_par[i, ], locgeo_optim_fn,
      gr = NULL,
      x = x, y = y, w = w,
      method = "L-BFGS-B",
      lower = c(0, 0, -max_speed, -max_speed),
      upper = c(pi, 2 * pi, max_speed, max_speed),
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

.estimate_locgeo <- function(x, y, x_new, kernel_fun, h, max_speed, grid_size, accuracy) {
  estim <- matrix(NA_real_, ncol=3, nrow=length(x_new))
  for (i in seq_along(x_new)) {
    t <- x_new[i]
    w <- kernel_fun((x-t)/h) / h
    w <- w / sum(w)
    if (!all(is.finite(w))) stop("Not all weights are finite.")
    res <- locgeo_optim(x, y, w, max_speed, grid_size, accuracy)
    estim[i,] <- Exp(res$p, t*res$v)
  }
  estim_a <- convert_e2a(estim)
  list(estim=estim, estim_a=estim_a)
}


#' Local geodesic regression on the sphere with cross validation.
#'
#' Uses leave one out cross validation to find a suitable bandwidth
#'
#' @param bw number of bandwidths to check
#' @param accuracy double in [0,1]. trade optimization accuracy for
#'   computational speed
#' @export
estimate_locgeo <- function(
  x, y, x_new,
  adapt=c("loocv", "none"),
  bw=7, kernel = "epanechnikov",
  max_speed=10, grid_size = 2, accuracy=0.25
) {
  kernel_fun <- get_kernel_fun(kernel)
  adapt <- match.arg(adapt)
  if (adapt == "loocv") {
    n <- length(x)
    hs <- (2/n)^seq(1, 0, len=bw)
    dists <- sapply(hs, function(h) {
      v <- sapply(seq_along(x), function(j) {
        res <- .estimate_locgeo(x[-j], y[-j,], x[j],
                                kernel_fun, h, max_speed, grid_size, accuracy)
        dist(res$estim, y[j, ])
      })
      mean(v)
    })
    h <- hs[which.min(dists)]
  } else {
    h <- bw
  }
  c(.estimate_locgeo(x, y, x_new, kernel_fun, h, max_speed, grid_size, accuracy), list(h = h))
}
