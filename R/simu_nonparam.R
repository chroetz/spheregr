#' @export
nonparam_methods <- c("locfre", "trifre", "locgeo", "trigeo")

#' @export
nonparam_colors <- list(
  locfre = rgb(1, 0, 0),
  trifre = rgb(0, 1, 0),
  locgeo = rgb(0, 0, 1),
  trigeo = rgb(1, 0, 1)
)

spiral <- function(n,
                   theta_min,
                   theta_max = theta_min,
                   phi_start = 0.5,
                   circles = 1) {
  phi <- seq(phi_start, phi_start + 2 * pi * circles, len = n) %% (2 * pi)
  theta <- seq(theta_min, theta_max, len = n)
  m_a <- cbind(theta, phi)
  m_a
}


sample_spiral <- function(n, n_new, sd, ...) {
  m_a <- spiral(n, ...)
  m_new_a <- spiral(n_new, ...)

  x_new <- seq(0, 1, len = n_new)
  x <- seq(0, 1, len = n)
  m <- convert_a2e(m_a)
  contr <- max(0, min(1, sd / sqrt((pi ^ 2 - 4) / 2)))
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
    m_new_a = m_new_a,
    sd = sd,
    contr = contr
  )
}



nonparam_reg <- function(method = nonparam_methods, ...) {
  method <- match.arg(method)
  switch(method,
    locfre = estimate_locfre(...),
    trifre = estimate_trifre(...),
    locgeo = estimate_locgeo(...),
    trigeo = estimate_trigeo(...)
  )
}

nonparam_run_one <- function(osamp, ometh) {
  data <- do.call(sample_spiral, osamp)
  res <- list()
  opt_data <- list(x = data$x,
                   y = data$y,
                   x_new = data$x_new)
  for (meth in linreg_methods) {
    res[[meth]] <-
      do.call(nonparam_reg, c(ometh[[meth]], opt_data, list(method = meth)))
  }
  c(data, res)
}

create_nonparam_opts <- function(reps, n, sd, curve = c("simple", "spiral")) {
  o <- list(samp = list(),
            simu = list(),
            meth = list())

  o$simu$reps <- reps

  o$samp$n <- n
  o$samp$n_new <- 100
  o$samp$noise <- "contracted"
  o$samp$sd <- sd

  curve <- match.arg(curve)
  switch(curve,
    simple = {
      o$samp$theta_min <- pi/4
      o$samp$theta_max <- pi/4
      o$samp$phi_start <- 0.5
      o$samp$circles <- 1
    },
    spiral = {
      o$samp$theta_min <- pi/8
      o$samp$theta_max <- pi/8*7
      o$samp$phi_start <- 0.5
      o$samp$circles <- 1.5
    }
  )

  o$meth[nonparam_methods] <- list(list())

  o$meth$locfre$loocv <- TRUE
  o$meth$locfre$kernel <- "epanechnikov"
  o$meth$locfre$bw <- 7
  o$meth$locfre$restarts <- 2

  o$meth$trifre$adapt <- "loocv"
  o$meth$trifre$num_basis <- 20
  o$meth$trifre$periodize <- FALSE
  o$meth$trifre$restarts <-  2

  o$meth$locgeo$adapt <- "loocv"
  o$meth$locgeo$bw <- 7
  o$meth$locgeo$kernel <- "epanechnikov"
  o$meth$locgeo$max_speed <- 10
  o$meth$locgeo$restarts <- NULL
  o$meth$locgeo$accuracy <- 0.3

  o$meth$trigeo$adapt <- "none"
  o$meth$trigeo$num_basis <- 3
  o$meth$trigeo$periodize <- FALSE
  o$meth$trigeo$max_speed <- 5
  o$meth$trigeo$restarts <- NULL
  o$meth$trigeo$accuracy <- 0.25

  o
}


nonparam_simulate <- function(opts) {
  lapply(opts, function(opt) {
    replicate(opt$simu$reps, nonparam_run_one(opt$samp, opt$meth), simplify = FALSE)
  })
}
