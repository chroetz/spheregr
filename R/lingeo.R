#' @param p 1x3
#' @param v 1x3
#' @param x nx1
#' @param y nx3
energy <- function(p, v, x, y) {
  z <- x %*% v # nx3
  n <- nrow(z)
  norm_z <- sqrt(.rowSums(z^2, n, 3L, FALSE))
  co <- tcrossprod(y, p) * cos(norm_z) # nx1
  si <- .rowSums(y * z, n, 3L, FALSE) / norm_z * sin(norm_z) # nx1s
  s <- clamp1(co + si)
  .Internal(mean(acos(s)^2)) / 2
  # should be same as
  # res2 <- mean(dist_sphere_R3(Exp(p, z), y)^2) / 2
}

optim_fn <- function(par, x, y) {
  p <- convert_a2e_1(par[1:2])
  v <- get_v_rcpp(p, par[3:4])
  dim(p) <- c(1L, 3L)
  dim(v) <- c(1L, 3L)
  energy(p, v, x, y)
}

exec_optim <- function(x, y, max_speed, grid_size, accuracy=0.5) {
  initial_parameters <- get_initial_parameters(grid_size, max_speed)
  res_lst <- list()
  for (i in seq_len(nrow(initial_parameters))) {
    res_lst[[i]] <- stats::optim(
      initial_parameters[i, ], optim_fn,
      gr = NULL,
      x = x, y = y,
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

#' Use geodesic regression to find the best fitting geodesic
#'
#' @param x covariate, a vector of length n with values in [0, 1]
#' @param y observations on the sphere (in R3), a nx3 matrix
#' @param x_new double vector, new values of the covariate where to evaluate the
#'   estimated function
#' @param max_speed a nonnegative scalar double, a bound on the maximum speed of
#'   the geodesic
estimate_lingeo <- function(x, y, x_new, max_speed, grid_size = 2) {
  res <- exec_optim(x, y, max_speed, grid_size)
  estim <- Exp(res$p, x_new %*% res$v)
  estim_a <- convert_e2a(estim)
  c(list(estim=estim, estim_a=estim_a), res)
}
