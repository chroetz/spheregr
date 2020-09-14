frechet_objective <- function(par, y_angle, yo, B_inv, X, X_eval_t) {
  yq <- acos(cos(y_angle[,1])*cos(par[1]) +
             sin(y_angle[,1])*sin(par[1])*cos(par[2]-y_angle[,2]))
  beta_hat_q <- B_inv %*% X %*% (yq^2 - yo^2)
  X_eval_t %*% beta_hat_q
}

prof <- profvis::profvis({
  res <- replicate(100, {
    speed_max <- 10
    n <- 30
    data <- sample_data(n=n, sd=0.3, speed_max=speed_max)
    x <- data$x
    y <- data$y
    y_angle <- R32angle(y)
    x_true <- x
    y_true <- data$y_true
    y_true_angle <- R32angle(y_true)
    p_true <- data$p
    v_true <- data$v

    initial_parameters <-
      expand.grid(
        alpha = (0:4) * pi/5 + pi/10,
        phi = (0:4) * 2*pi/5 + 2*pi/10
      ) %>%
      as.matrix()

    estim_angle <- matrix(nrow=n, ncol=2)

    X <- rbind(1, x)
    B_inv <- solve(tcrossprod(X))
    yo <- dist_angle(y_angle, matrix(c(0,0), nrow=1))

    for (j in 1:n) {
      X_eval_t = cbind(1, x_true[j])
      res_lst <- list()
      for (i in seq_len(nrow(initial_parameters))) {
        res_lst[[i]] <- optim(
          initial_parameters[i, ], frechet_objective, gr = NULL,
          X = X, y_angle = y_angle, X_eval_t=X_eval_t, yo=yo, B_inv=B_inv,
          method = "L-BFGS-B",
          lower = c(0, 0),
          upper = c(pi, 2*pi))
      }
      values <- sapply(res_lst, function(x) x$value)
      idx <- which.min(values)
      res <- res_lst[[idx]]
      estim_angle[j, ] <- res$par
    }
  })
})
