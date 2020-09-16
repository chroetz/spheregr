load("sim_rand")

# mise
sim_extract_mise(sim_rand)

# speed
sim_extract_speed_se(sim_rand)





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



load("sim_line")

# bias variance plots

s <- sim_line[[11]]
estims_a <- array(NA, dim=c(
  samples = length(s[[1]]$x_new),
  dims = 2,
  reps = length(s))
)
plot(NA,
     xlim=c(0, 2*pi), ylim=c(0, pi),
     xlab="phi", ylab="alpha")
sphere_grid()
striped_lines(s[[1]]$m_a[,2:1], col="gray", lwd=3)
for (meth in linreg_methods) {
  for (i in seq_along(s)) {
    estims_a[,,i] <- s[[i]][[meth]]$estim_a
  }
  mean_a_and_var <- t(apply(estims_a, 1, function(x) frechet_mean(t(x))))
  mean_a <- mean_a_and_var[,1:2]
  geo_sd <- sqrt(mean_a_and_var[,3])

#  plot(s[[1]]$x_new, geo_sd)
  striped_lines(mean_a[,2:1], col=method_cols[[meth]], lwd=3)
}

# for (i in 1:10) {
#   striped_lines(geodesic_a[,,i], col=rgb(0,1,0,0.5), lwd=1)
# }

plot(s[[1]]$x_new, s[[1]]$m_a[,1])
plot(s[[1]]$x_new, s[[1]]$m_a[,2])

