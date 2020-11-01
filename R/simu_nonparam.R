#' @export
nonparam_methods <- c("locfre", "trifre", "locgeo", "trigeo")

#' @export
nonparam_colors <- list(
  locfre = rgb(1,0,0),
  trifre = rgb(0,1,0),
  locgeo = rgb(0,0,1),
  trigeo = rgb(1,0,1)
)

ker_epan <- function(x) 3/4*(1-x^2)*(abs(x)<1)

spiral <- function(n, theta_min, theta_max=theta_min, phi_start = 0.5, circles=1) {
  phi <- seq(phi_start, phi_start+2*pi*circles, len=n) %% (2*pi)
  theta <- seq(theta_min, theta_max, len=n)
  m_a <- cbind(theta, phi)
  m_a
}


sample_spiral <- function(n, n_new, contr, ...) {
  n <- 40
  n_new <- 200
  contr <- 0.3

  m_a <- spiral(n, ...)
  m_new_a <- spiral(n_new, ...)

  x_new <- seq(0, 1, len=n_new)
  x <- seq(0, 1, len=n)
  m <- convert_a2e(m_a)
  sd <- sqrt((pi^2-4)/2)*contr
  y <- add_noise_contract(m, contr)
  y_a <- convert_e2a(y)

  list(
    n = n,
    n_new = n_new,
    x = x,
    x_new = x_new,
    y = y,
    y_a = y_a,
    m = m,
    m_a = m_a,
    sd = sd,
    contr = contr)
}



nonparam_reg <- function(x, y, x_new,
                   method = nonparam_methods,
                   ...) {
  method <- match.arg(method)
  switch(method,
   locfre = estimate_locfre(x, y, x_new, ...),
   trifre = estimate_trifre(x, y, x_new, ...),
   locgeo = estimate_locgeo(x, y, x_new, ...),
   trigeo = estimate_trigeo(x, y, x_new, ...)
  )
}

nonparam_run_one <- function(opt) {
  opt_spiral <- opt[]
  opt_reg <- opt[]
  data <- do.call(sample_spiral, opt_spiral)
  res <- list()
  opt_data <- list(x = data$x, y = data$y, x_new = data$x_new)
  for (meth in linreg_methods) {
    res[[meth]] <- do.call(nonparam_reg, c(opt_reg, opt_data, list(method = meth)))
  }
  c(data, res)
}


# TODO
create_nonparam_opts <- function(
  reps,
  x_new = seq(0, 1, len=300),
  max_speed = 10,
  n = c(10, 100)
) {
  ot <- expand.grid(
    n = n,
    sd = sd)
  linreg_simul_opts <- lapply(1:nrow(ot), function(i)
    list(
      n = ot[i,]$n,
      sd=ot[i,]$sd,
      speed=ot[i,]$speed,
      x_new=x_new,
      max_speed=max_speed,
      reps=reps))
}


nonparam_simulate <- function(opts) {
  lapply(opts, function(opt) {
    replicate(opt$reps, run_one(opt), simplify=FALSE)
  })
}
