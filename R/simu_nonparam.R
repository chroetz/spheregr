#' @export
nonparam_methods <- c("locfre", "trifre", "locgeo", "trigeo")

#' @export
nonparam_colors <- c(
  locfre = "#FF0000",
  trifre = "#00FF00",
  locgeo = "#0000FF",
  trigeo = "#FF00FF"
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


sample_spiral <- function(n, n_new, sd, noise, ...) {
  m_a <- spiral(n, ...)
  m_new_a <- spiral(n_new, ...)

  x_new <- seq(0, 1, len = n_new)
  x <- seq(0, 1, len = n)
  m <- convert_a2e(m_a)
  contr <- max(0, min(1, sd / sqrt((pi ^ 2 - 4) / 2)))
  y <- switch(noise,
              contracted = add_noise_contract(m, contr),
              normal = add_noise_normal(m, sd)
             )
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

nonparam_run_one <- function(osamp, ometh, verbosity=1) {
  data <- do.call(sample_spiral, osamp)
  res <- list()
  opt_data <- data[c("x", "y", "x_new")]
  for (meth in names(ometh)) {
    if (verbosity > 2) {
      cat("\t\t", meth, " ", sep="")
      pt <- proc.time()
    }
    res[[meth]] <-
      do.call(nonparam_reg, c(ometh[[meth]], opt_data, list(method = meth)))
    if (verbosity > 2) {
      cat((proc.time() - pt)[3], "\n")
    }
  }
  list(data=data, predict=res)
}

create_nonparam_opt <- function(
  reps, n, sd, n_new, curve = c("simple", "spiral"),
  trigeo=3, others=TRUE) {
  o <- list(samp = list(),
            simu = list(),
            meth = list())

  o$simu$reps <- reps

  o$samp$n <- n
  o$samp$n_new <- n_new
  o$samp$noise <- "contracted"
  o$samp$sd <- sd

  curve <- match.arg(curve)
  o$simu$curve <- curve
  switch(curve,
    simple = {
      o$samp$theta_min <- pi/4
      o$samp$theta_max <- pi/4
      o$samp$phi_start <- 0.5
      o$samp$circles <- 1
      periodic <- TRUE
    },
    spiral = {
      o$samp$theta_min <- pi/8
      o$samp$theta_max <- pi/8*7
      o$samp$phi_start <- 0.5
      o$samp$circles <- 1.5
      periodic <- FALSE
    }
  )

  if (others) {
    o$meth$locfre$adapt <- "loocv"
    o$meth$locfre$kernel <- "epanechnikov"
    o$meth$locfre$bw <- 7
    o$meth$locfre$restarts <- 3

    o$meth$trifre$adapt <- "loocv"
    o$meth$trifre$num_basis <- 20
    o$meth$trifre$periodize <- !periodic
    o$meth$trifre$restarts <-  3

    o$meth$locgeo$adapt <- "loocv"
    o$meth$locgeo$bw <- 7
    o$meth$locgeo$kernel <- "epanechnikov"
    o$meth$locgeo$max_speed <- 5
    o$meth$locgeo$restarts <- NULL
    o$meth$locgeo$accuracy <- 0.25
  }
  if (trigeo) {
    o$meth$trigeo$adapt <- "none"
    o$meth$trigeo$num_basis <- trigeo
    o$meth$trigeo$periodize <- !periodic
    o$meth$trigeo$max_speed <- 5
    o$meth$trigeo$restarts <- NULL
    o$meth$trigeo$accuracy <- 0.25
  }

  o
}


nonparam_simulate <- function(opt_list, verbosity=1) {
  all_res <- list()
  for (i in seq_along(opt_list)) {
    opt <- opt_list[[i]]
    if (verbosity > 0) {
      cat("start opt", i,"/", length(opt_list),"\n")
      cat("\truns: ", opt$simu$reps, "\n")
      cat("\tn: ", opt$samp$n, "\n")
      pto <- proc.time()
    }
    res <- list()
    for (j in seq_len(opt$simu$reps)) {
      if (verbosity > 1) {
        cat("\tstart run", j,"/", opt$simu$reps, "\n")
        pt <- proc.time()
      }
      res[[j]] <- nonparam_run_one(opt$samp, opt$meth, verbosity)
      if (verbosity > 1) {
        cat("\t\ttotal:", (proc.time() - pt)[3], "\n")
        cat("\tend run", j,"/", opt$simu$reps,"\n")
      }
    }
    all_res[[i]] <- res
    if (verbosity > 0) {
      cat("\ttime:", (proc.time() - pto)[3], "\n")
      cat("end opt", i,"/", length(opt_list),"\n")
    }
  }
  all_res
}


nonparam_simulate_parallel <- function(opt_list) {
  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(cores)
  parallel::clusterEvalQ(cl, devtools::load_all(".")) # TODO!
  all_res <- lapply(opt_list, function(opt) {
    parallel::parLapply(cl, seq_len(opt$simu$reps), nonparam_run_one, osamp=opt$samp, ometh=opt$meth)
  })
  parallel::stopCluster(cl)
  all_res
}
