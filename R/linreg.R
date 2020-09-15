sample_data <- function(n, sd, speed_max) {
  x <- (1:n) / n
  p_angle <- runif(2) * c(pi, 2 * pi)
  p <- angle2R3(p_angle)
  speed <- runif(1, min = 0, max = speed_max)
  v <- norm_vec(t(pracma::nullspace(p) %*% rnorm(2))) * speed
  y_true <- Exp(p, x %*% v)
  y <- y_true
  for (i in 1:n) {
    noise_dir <- norm_vec(t(pracma::nullspace(y_true[i, , drop = FALSE]) %*% rnorm(2)))
    y[i, ] <- Exp(y_true[i, ], noise_dir * rnorm(1, sd = sd))
  }
  list(x = x, y = y, p = p, v = v, y_true = y_true, speed = speed)
}



linreg <- function(x, y, x_new,
                   method = c("frechet", "geodesic", "cosine"),
                   restarts = 2, max_speed = 10) {
  method <- match.arg(method)
  switch(method,
    frechet = estimate_frechet(x, y, x_new, restarts),
    geodesic = estimate_geodesic(x, y, x_new, restarts, max_speed = max_speed),
    cosine = estimate_cosine(x, y, x_new, restarts, max_speed = max_speed)
  )
}
