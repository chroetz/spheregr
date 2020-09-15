sample_data <- function(n, sd, x_eval=NULL, speed_bounds) {
  x <- (1:n) / n
  if (is.null(x_eval)) x_eval <- x
  p_a <- matrix(runif(2) * c(pi, 2 * pi), ncol=2)
  p <- convert_a2e(p_a)
  speed <- runif(1, min = speed_bounds[1], max = speed_bounds[2])
  v <- norm_vec(t(pracma::nullspace(p) %*% rnorm(2))) * speed
  y <- Exp(p, x %*% v)
  for (i in 1:n) {
    noise_dir <- norm_vec(t(pracma::nullspace(y[i, , drop = FALSE]) %*% rnorm(2)))
    y[i, ] <- Exp(y[i, ], noise_dir * rnorm(1, sd = sd))
  }
  m <- Exp(p, x_eval %*% v)
  list(
    x = x,
    x_eval = x_eval,
    y = y,
    y_a = convert_e2a(y),
    p = p,
    v = v,
    m = m,
    m_a = convert_e2a(m),
    speed = speed)
}



linreg <- function(x, y, x_new,
                   method = c("frechet", "geodesic", "cosine"),
                   speed_bounds = c(0, 10),
                   restarts = 2) {
  method <- match.arg(method)
  switch(method,
    frechet = estimate_frechet(x, y, x_new, restarts),
    geodesic = estimate_geodesic(x, y, x_new, speed_bounds, restarts),
    cosine = estimate_cosine(x, y, x_new, speed_bounds, restarts)
  )
}
