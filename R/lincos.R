#' Estimate lambda of the (1D) cosine regression.
#'
#' @param x double(n) in [0,1], covariates.
#' @param y nx3 matrix, observations on sphere.
#' @param q kx3 matrix, test points on sphere.
#' @return double(1), estimated lambda
cos_estim_speed <- function(x, y, q, max_speed) {
  cos_dist_yq <- y %*% t(q)
  objective <- function(par) {
    X <- rbind(cos(par * x), sin(par * x))
    B <- X %*% t(X)
    pars_ab <- solve(B, X) %*% cos_dist_yq
    mean((t(X) %*% pars_ab - cos_dist_yq)^2)
  }
  stats::optimize(objective, lower = 0, upper = max_speed)
}


cos_objective <- function(q_a, BinvXy, x_new, speed) {
  q <- convert_a2e_1(q_a)
  pars_ab <- BinvXy %*% q
  -(pars_ab[1] * cos(speed * x_new) + pars_ab[2] * sin(speed * x_new))
}


estimate_lincos <- function(x, y, x_new, max_speed, grid_size = 2) {
  # estimate speed
  q_a <- get_initial_parameters(grid_size, num_basis=0)
  q <- convert_a2e(q_a)
  res <- cos_estim_speed(x, y, q, max_speed)
  speed <- res$minimum

  # estimate m(t) via optim()
  X <- rbind(cos(speed * x), sin(speed * x))
  B <- X %*% t(X)
  BinvXy <- solve(B, X %*% y)
  initial_parameters <-  get_initial_parameters(grid_size, num_basis=0)
  estim_a <- matrix(nrow = length(x_new), ncol = 2)
  for (j in seq_along(x_new)) {
    res_lst <- list()
    for (i in seq_len(nrow(initial_parameters))) {
      res_lst[[i]] <- stats::optim(
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
