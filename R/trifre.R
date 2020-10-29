trifre_objective <- function(par, y_a, yo, B) {
  yq <- acos(cos(y_a[, 1]) * cos(par[1]) +
    sin(y_a[, 1]) * sin(par[1]) * cos(par[2] - y_a[, 2]))
  m_hat <- B %*% (yq^2 - yo^2)
}

eval_trig_base <- function(x, n_basis) {
  k <- seq_len(n_basis/2)
  phi <- matrix(NA, nrow=2*length(k)+1, ncol=length(x))
  phi[1, ] <- 1
  phi[2*k, ] <- sqrt(2)*cospi(2 * k %*% t(x))
  phi[2*k+1, ] <- sqrt(2)*sinpi(2 * k %*% t(x))
  phi[1:n_basis, ]
}

#' @export
estimate_trifre <- function(x, y, x_new, n_basis, restarts = 2) {
  N <- restarts
  initial_parameters <-
    expand.grid(
      alpha = (0:(N - 1)) * pi / N + pi / N / 2,
      phi = (0:(N - 1)) * 2 * pi / N + 2 * pi / N / 2
    ) %>%
    as.matrix()

  estim_a <- matrix(nrow = length(x_new), ncol = 2)

  y_a <- convert_e2a(y)
  X <- rbind(1, x)
  yo <- dist_a(y_a, matrix(c(0, 0), nrow = 1))

  for (j in seq_along(x_new)) {
    res_lst <- list()
    phi <- eval_trig_base(x, n_basis)
    phi_new <- eval_trig_base(x_new[j], n_basis)
    B <- 1/n * crossprod(phi_new, phi)
    for (i in seq_len(nrow(initial_parameters))) {
      res_lst[[i]] <- stats::optim(
        initial_parameters[i, ], trifre_objective,
        gr = NULL,
        y_a = y_a, B = B, yo = yo,
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
