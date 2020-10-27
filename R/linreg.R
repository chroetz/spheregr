#' @export
linreg_methods <- c("frechet", "geodesic", "cosine")

#' @export
method_colors <- list(
  frechet = rgb(1,0,0),
  geodesic = rgb(0,1,0),
  cosine = rgb(0,0,1))


#' Sample regression data from a geodesic model with normal noise.
#'
#' @param n positive integer. number of observations to sample.
#' @param sd positive double. standard deviation of normal noise.
#' @param x_new double vector, new values of the covariate where to evaluate the
#'   estimated function
#' @param speed_bounds double vector of length 1 or 2
#' @return a list with observations at x in y and true function values at x_new
#'   in m
sample_data <- function(n, sd, x_new=NULL, speed_bounds=NULL, p=NULL, v=NULL) {
  x <- (1:n) / n
  if (is.null(x_new)) x_new <- x
  if (!is.null(speed_bounds)) {
    speed_target <- if (length(speed_bounds) == 2) {
      runif(1, min = speed_bounds[1], max = speed_bounds[2])
    } else {
      speed_bounds
    }
  } else {
    speed_target <- NULL
  }
  if (is.null(p) || is.null(v)) {
    p_a <- matrix(runif(2) * c(pi, 2 * pi), ncol=2)
    p <- convert_a2e(p_a)
    v <- norm_vec(t(pracma::nullspace(p) %*% rnorm(2))) * speed_target
  } else {
    p_a <- convert_e2a(p)
    if (!is.null(speed_target)) v <- v * speed_target / sqrt(sum(v^2))
  }
  speed <- sqrt(sum(v^2))
  m <- Exp(p, x_new %*% v)
  y <- add_noise_normal(Exp(p, x %*% v), sd)

  list(
    x = x,
    x_new = x_new,
    y = y,
    y_a = convert_e2a(y),
    p = p,
    p_a = p_a,
    v = v,
    m = m,
    m_a = convert_e2a(m),
    speed = speed,
    n = n,
    sd = sd,
    speed_bounds = speed_bounds)
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


run_one <- function(opt) {
  data <- sample_data(opt$n, opt$sd, x_new=opt$x_new, speed_bounds=opt$speed, v=opt$v, p=opt$p)
  res <- list()
  for (meth in linreg_methods) {
    res[[meth]] <- linreg(data$x, data$y, data$x_new, method=meth, speed_bounds=opt$speed_bounds)
  }
  c(data, res)
}

create_linreg_opts <- function(
  reps,
  x_new = seq(0, 1, len=300),
  speed_bounds = c(0, 10),
  n = c(10, 100),
  sd = c(0.1, 1),
  speed = c(1, pi, 8),
  p = NULL,
  v = NULL
) {
  ot <- expand.grid(
    n = n,
    sd = sd,
    speed = speed)
  linreg_simul_opts <- lapply(1:nrow(ot), function(i)
    list(
      n = ot[i,]$n,
      sd=ot[i,]$sd,
      speed=ot[i,]$speed,
      x_new=x_new,
      speed_bounds=speed_bounds,
      reps=reps,
      p = p,
      v = v))
}

simulate <- function(opts) {
  lapply(opts, function(opt) {
    replicate(opt$reps, run_one(opt), simplify=FALSE)
  })
}
