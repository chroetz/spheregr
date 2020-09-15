frechet_objective <- function(par, y_angle, yo, B_inv, X, X_eval_t) {
  yq <- acos(cos(y_angle[, 1]) * cos(par[1]) +
    sin(y_angle[, 1]) * sin(par[1]) * cos(par[2] - y_angle[, 2]))
  beta_hat_q <- B_inv %*% X %*% (yq^2 - yo^2)
  X_eval_t %*% beta_hat_q
}

estimate_frechet <- function(x, y, x_new, restarts = 2) {
  N <- restarts
  initial_parameters <-
    expand.grid(
      alpha = (0:(N - 1)) * pi / N + pi / N / 2,
      phi = (0:(N - 1)) * 2 * pi / N + 2 * pi / N / 2
    ) %>%
    as.matrix()

  estim_angle <- matrix(nrow = length(x_new), ncol = 2)

  y_angle <- R32angle(y)
  X <- rbind(1, x)
  B_inv <- solve(tcrossprod(X))
  yo <- dist_angle(y_angle, matrix(c(0, 0), nrow = 1))

  for (j in seq_along(x_new)) {
    X_eval_t <- cbind(1, x_new[j])
    res_lst <- list()
    for (i in seq_len(nrow(initial_parameters))) {
      res_lst[[i]] <- optim(
        initial_parameters[i, ], frechet_objective,
        gr = NULL,
        X = X, y_angle = y_angle, X_eval_t = X_eval_t, yo = yo, B_inv = B_inv,
        method = "L-BFGS-B",
        lower = c(0, 0),
        upper = c(pi, 2 * pi)
      )
    }
    values <- sapply(res_lst, function(x) x$value)
    idx <- which.min(values)
    res <- res_lst[[idx]]
    estim_angle[j, ] <- res$par
  }

  angle2R3(estim_angle)
}
