#' Estimate lambda of the (1D) cosine regression.
#'
#' @param x double(n) in [0,1], covariates.
#' @param y nx3 matrix, observations on sphere.
#' @param q kx3 matrix, test points on sphere.
#' @return double(1), estimated lambda
cos_estim_speed <- function(x, y, q, speed_bounds) {
  cos_dist_yq <- y %*% t(q)
  objective <- function(par) {
    X <- rbind(cos(par * x), sin(par * x))
    B <- X %*% t(X)
    pars_ab <- solve(B, X) %*% cos_dist_yq
    mean((t(X) %*% pars_ab - cos_dist_yq)^2)
  }
  optimize(objective, speed_bounds)
}


cos_objective <- function(q_a, BinvXy, x_new, speed) {
  q <- convert_a2e_1(q_a)
  pars_ab <- BinvXy %*% q
  -(pars_ab[1] * cos(speed * x_new) + pars_ab[2] * sin(speed * x_new))
}


estimate_cosine <- function(x, y, x_new, speed_bounds, restarts = 2) {
  # estimate speed
  N <- 2
  M <- 2
  q_a <- expand.grid(
    alpha = (0:(N - 1)) * pi / N + pi / N / 2,
    phi = (0:(M - 1)) * 2 * pi / M + 2 * pi / M / 2
  ) %>%
    as.matrix()
  q <- convert_a2e(q_a)
  res <- cos_estim_speed(x, y, q, speed_bounds)
  speed <- res$minimum

  # estimate m(t) via optim()
  X <- rbind(cos(speed * x), sin(speed * x))
  B <- X %*% t(X)
  BinvXy <- solve(B, X %*% y)
  N <- restarts
  initial_parameters <- expand.grid(
    alpha = (0:(N - 1)) * pi / N + pi / N / 2,
    phi = (0:(N - 1)) * 2 * pi / N + 2 * pi / N / 2
  ) %>%
    as.matrix()
  estim_a <- matrix(nrow = length(x_new), ncol = 2)
  for (j in seq_along(x_new)) {
    res_lst <- list()
    for (i in seq_len(nrow(initial_parameters))) {
      res_lst[[i]] <- optim(
        initial_parameters[i, ], cos_objective,
        gr = NULL,
        BinvXy = BinvXy, x_new = x_new[j], speed = speed,
        method = "L-BFGS-B",
        lower = c(0, 0),
        upper = c(pi, 2 * pi)
      )
    }
    values <- sapply(res_lst, function(x) x$value)
    idx <- which.min(values)
    res <- res_lst[[idx]]
    estim_a[j, ] <- res$par
  }

  estim <- convert_a2e(estim_a)
  list(estim=estim, estim_a=estim_a, speed=speed)
}
