locfre_objective <- function(par, y_a, yo, w, X, X_eval_t) {
  B_inv <- solve(X %*% diag(w) %*% t(X))
  yq <- acos(cos(y_a[, 1]) * cos(par[1]) +
    sin(y_a[, 1]) * sin(par[1]) * cos(par[2] - y_a[, 2]))
  beta_hat_q <- B_inv %*% X %*% diag(w) %*% (yq^2 - yo^2)
  X_eval_t %*% beta_hat_q
}

.estimate_locfre <- function(x, y, x_new, kernel, h, restarts) {
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
    t <- x_new[j]
    X_eval_t <- cbind(1, t)
    w <- kernel((x-t)/h) / h
    w <- w / sum(w)
    if (!all(is.finite(w))) stop("Not all weights are finite.")
    res_lst <- list()
    for (i in seq_len(nrow(initial_parameters))) {
      res_lst[[i]] <- stats::optim(
        initial_parameters[i, ], locfre_objective,
        gr = NULL,
        X = X, y_a = y_a, X_eval_t = X_eval_t, yo = yo, w = w,
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

#' Local linear Frechet regression on the sphere with cross validation.
#'
#' Uses leave one out cross validation to find a suitable bandwidth
#'
#' @param bw the bandwidth $h$ or number of bandwidths to check (for adapt="loocv")
#' @adapt either "loocv" or NA
#' @export
estimate_locfre <- function(
  x, y, x_new,
  adapt=c("loocv", "none"), bw=7, kernel="epanechnikov", restarts = 2
) {
  kernel <- get_kernel_fun(kernel)
  adapt <- match.arg(adapt)
  if (adapt == "loocv") {
    n <- length(x)
    hs <- (3/n)^seq(1, 0, len=bw)
    dists <- sapply(hs, function(h) {
      v <- sapply(seq_along(x), function(j) {
        res <- .estimate_locfre(x[-j], y[-j,], x[j], h=h, ...)
        dist(res$estim, y[j, ])
      })
      mean(v)
    })
    h <- hs[which.min(dists)]
  } else {
    h <- bw
  }
  c(.estimate_locfre(x, y, x_new, h=h, ...), list(h = h))
}
