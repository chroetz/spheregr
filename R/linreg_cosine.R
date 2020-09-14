#' Estimate lambda of the (1D) cosine regression.
#'
#' @param x double(n) in [0,1], covariates.
#' @param y nx3 matrix, observations on sphere.
#' @param q kx3 matrix, test points on sphere.
#' @return double(1), estimated lambda
cos_estim_speed <- function(x, y, q, max_speed=10) {
  cos_dist_yq <- y %*% t(q)
  objective <- function(par) {
    X <- rbind(cos(par*x), sin(par*x))
    B <- X %*% t(X)
    pars_ab <- solve(B, X) %*% cos_dist_yq
    mean((t(X) %*% pars_ab - cos_dist_yq)^2)
  }
  optimize(objective, c(0, max_speed))
}


cos_objective <- function(q_angle, BinvXy, x_new, speed) {
  q <- angle2R3_1(q_angle)
  pars_ab <- BinvXy %*% q
  -(pars_ab[1]*cos(speed*x_new) + pars_ab[2]*sin(speed*x_new))
}


estimate_cosine <- function(x, y, x_new, restarts=2, max_speed=10) {
  # estimate speed
  N <- 2
  M <- 2
  q_angle <- expand.grid(
    alpha = (0:(N-1)) * pi/N + pi/N/2,
    phi = (0:(M-1)) * 2*pi/M + 2*pi/M/2
  ) %>%
    as.matrix()
  q <- angle2R3(q_angle)
  res <- cos_estim_speed(x, y, q, max_speed=max_speed)
  speed <- res$minimum

  # estimate m(t) via optim()
  X <- rbind(cos(speed*x), sin(speed*x))
  B <- X %*% t(X)
  BinvXy <- solve(B, X %*% y)
  N <- restarts
  initial_parameters <- expand.grid(
      alpha = (0:(N-1)) * pi/N + pi/N/2,
      phi = (0:(N-1)) * 2*pi/N + 2*pi/N/2
    ) %>%
    as.matrix()
  estim_angle <- matrix(nrow=length(x_new), ncol=2)
  for (j in seq_along(x_new)) {
    res_lst <- list()
    for (i in seq_len(nrow(initial_parameters))) {
      res_lst[[i]] <- optim(
        initial_parameters[i, ], cos_objective, gr = NULL,
        BinvXy = BinvXy, x_new = x_new[j], speed = speed,
        method = "L-BFGS-B",
        lower = c(0, 0),
        upper = c(pi, 2*pi))
    }
    values <- sapply(res_lst, function(x) x$value)
    idx <- which.min(values)
    res <- res_lst[[idx]]
    estim_angle[j, ] <- res$par
  }

  angle2R3(estim_angle)
}

