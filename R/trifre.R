trifre_objective <- function(par, y_a, yo, B) {
  yq <- acos(cos(y_a[, 1]) * cos(par[1]) +
    sin(y_a[, 1]) * sin(par[1]) * cos(par[2] - y_a[, 2]))
  m_hat <- B %*% (yq^2 - yo^2)
}


.estimate_trifre <- function(x, y, x_new, num_basis, periodize, grid_size) {
  initial_parameters <- get_initial_parameters(grid_size, num_basis=0)
  estim_a <- matrix(nrow = length(x_new), ncol = 2)
  n <- length(x)

  if (periodize) {
    x <- c(x, rev(2-x))/2
    y <- rbind(y, y[n:1,])
    x_new <- x_new/2
  }

  y_a <- convert_e2a(y)
  X <- rbind(1, x)
  yo <- dist_a(y_a, matrix(c(0, 0), nrow = 1))

  for (j in seq_along(x_new)) {
    res_lst <- list()
    phi <- eval_trig_base(x, num_basis)
    phi_new <- eval_trig_base(x_new[j], num_basis)
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

#' @export
estimate_trifre <- function(x, y, x_new, adapt=c("loocv", "none"), num_basis=20,
                            periodize=FALSE, grid_size = 2) {
  adapt <- match.arg(adapt)
  if (adapt == "loocv") {
    n <- length(x)
    nums_basis <- unique(round(1 + 2^(seq(0, log2(n/2), len=num_basis))))
    dists <- sapply(nums_basis, function(num_b) {
      v <- sapply(seq_along(x), function(j) {
        res <- .estimate_trifre(
          x[-j], y[-j,], x[j],
          num_b, periodize, grid_size)
        dist(res$estim, y[j, ])
      })
      mean(v)
    })
    num_b <- nums_basis[which.min(dists)]
  } else {
    num_b <- num_basis
  }
  c(.estimate_trifre(x, y, x_new, num_b, periodize, grid_size), list(num_basis = num_b))
}
