frechet_objective <- function(par, y_a, yo, B_inv, X, X_eval_t) {
  yq <- acos(cos(y_a[, 1]) * cos(par[1]) +
    sin(y_a[, 1]) * sin(par[1]) * cos(par[2] - y_a[, 2]))
  beta_hat_q <- B_inv %*% X %*% (yq^2 - yo^2)
  X_eval_t %*% beta_hat_q
}

estimate_linfre <- function(x, y, x_new, grid_size = 2) {
  initial_parameters <- get_initial_parameters(grid_size, num_basis=0)
  estim_a <- matrix(nrow = length(x_new), ncol = 2)

  y_a <- convert_e2a(y)
  X <- rbind(1, x)
  B_inv <- solve(tcrossprod(X))
  yo <- dist_a(y_a, matrix(c(0, 0), nrow = 1))

  for (j in seq_along(x_new)) {
    X_eval_t <- cbind(1, x_new[j])
    res_lst <- list()
    for (i in seq_len(nrow(initial_parameters))) {
      res_lst[[i]] <- stats::optim(
        initial_parameters[i, ], frechet_objective,
        gr = NULL,
        X = X, y_a = y_a, X_eval_t = X_eval_t, yo = yo, B_inv = B_inv,
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
  list(estim=estim, estim_a=estim_a)
}
