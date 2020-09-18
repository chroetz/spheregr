
frechet_mean_objective <- function(par, y_a) {
  yq <- acos(cos(y_a[, 1]) * cos(par[1]) +
               sin(y_a[, 1]) * sin(par[1]) * cos(par[2] - y_a[, 2]))
  mean(yq^2)
}

frechet_mean <- function(y_a, restarts=2) {
  N <- restarts
  initial_parameters <-
    expand.grid(
      alpha = (0:(N - 1)) * pi / N + pi / N / 2,
      phi = (0:(N - 1)) * 2 * pi / N + 2 * pi / N / 2
    ) %>%
    as.matrix()
  res_lst <- list()
  for (i in seq_len(nrow(initial_parameters))) {
    res_lst[[i]] <- stats::optim(
      initial_parameters[i, ], frechet_mean_objective,
      gr = NULL,
      y_a = y_a,
      method = "L-BFGS-B",
      lower = c(0, 0),
      upper = c(pi, 2 * pi)
    )
  }
  values <- sapply(res_lst, function(x) x$value)
  idx <- which.min(values)
  res <- res_lst[[idx]]
  c(alpha=res$par[1], phi=res$par[2], var=res$value)
}
