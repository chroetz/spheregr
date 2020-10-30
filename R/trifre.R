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
estimate_trifre <- function(x, y, x_new, n_basis, periodize=FALSE, restarts = 2) {
  N <- restarts
  initial_parameters <-
    expand.grid(
      alpha = (0:(N - 1)) * pi / N + pi / N / 2,
      phi = (0:(N - 1)) * 2 * pi / N + 2 * pi / N / 2
    ) %>%
    as.matrix()

  estim_a <- matrix(nrow = length(x_new), ncol = 2)

  if (periodize) {
    n <- length(x)
    x <- c(x, rev(2-x))/2
    y <- rbind(y, y[n:1,])
    x_new <- x_new/2
  }

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


estimate_trifre_loocv <- function(x, y, x_new, n_n_basis=20, ...) {
  n <- length(x)
  ns_basis <- unique(round(1 + 2^(seq(0, log2(n/2), len=n_n_basis))))
  dists <- sapply(ns_basis, function(n_basis) {
    v <- sapply(seq_along(x), function(j) {
      res <- estimate_trifre(x[-j], y[-j,], x[j], n_basis=n_basis, ...)
      dist(res$estim, y[j, ])
    })
    mean(v)
  })
  n_basis <- ns_basis[which.min(dists)]
  c(estimate_trifre(x, y, x_new, n_basis=n_basis, ...), list(n_basis = n_basis))
}
