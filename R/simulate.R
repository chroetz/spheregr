#' @export
methods <- c("linfre", "lincos", "lingeo", "locfre", "trifre", "locgeo", "trigeo")


#' @export
nonparam_colors <- c(
  linfre = "#FFFF00",
  lincos = "#00FFFF",
  lingeo = "#FF7FFF",
  locfre = "#FF0000",
  trifre = "#00FF00",
  locgeo = "#0000FF",
  trigeo = "#FF00FF"
)


spiral_fun <- function(theta_min, theta_max = theta_min, phi_start = 0.5, circles = 1) {
  function(x) {
    phi <- (phi_start + x * 2 * pi * circles) %% (2 * pi)
    theta <- theta_min + x * (theta_max - theta_max)
    cbind(theta, phi)
  }
}


geodesic_fun <- function(speed_bounds=NULL, p=NULL, v=NULL) {
  speed_target <-
    if (length(speed_bounds) == 2) {
      runif(1, min = speed_bounds[1], max = speed_bounds[2])
    } else {
      speed_bounds
    }
  if (is.null(p) || is.null(v)) {
    p_a <- matrix(runif(2) * c(pi, 2 * pi), ncol=2)
    p <- convert_a2e(p_a)
    v <- norm_vec(t(pracma::nullspace(p) %*% rnorm(2))) * speed_target
  } else {
    if (!is.null(speed_target)) v <- v * speed_target / sqrt(sum(v^2))
  }
  function(x) Exp(p, x %*% v)
}



sample_regression_data <- function(n,
                                   n_new,
                                   sd,
                                   noise = c("contracted", "normal"),
                                   curve = c("geodesic", "spiral"),
                                   ...) {
  x_new <- seq(0, 1, len = n_new)
  x <- seq(0, 1, len = n)
  curve <- match.arg(curve)
  m_fun <- switch(curve,
                  geodesic = geodesic_fun(...),
                  spiral = spiral_fun(...))
  m_a <- m_fun(x, ...)
  m_new_a <- m_fun(x_new, ...)
  m <- convert_a2e(m_a)
  noise <- match.arg(noise)
  y <- switch(noise,
              contracted = add_noise_contract(m, sd = sd),
              normal = add_noise_normal(m, sd = sd))
  list(
    x = x,
    x_new = x_new,
    y = y,
    y_a = convert_e2a(y),
    m = m,
    m_a = m_a,
    m_new = convert_a2e(m_new_a),
    m_new_a = m_new_a
  )
}


estimate <- function(method = methods, ...) {
  method <- match.arg(method)
  switch(method,
         linfre = estimate_linfre(...),
         lincos = estimate_lincos(...),
         lingeo = estimate_lingeo(...),
         locfre = estimate_locfre(...),
         trifre = estimate_trifre(...),
         locgeo = estimate_locgeo(...),
         trigeo = estimate_trigeo(...)
  )
}


simu_run_one <- function(osamp, ometh, verbosity=1) {
  data <- do.call(sample_regression_data, osamp)
  res <- list()
  opt_data <- data[c("x", "y", "x_new")]
  for (meth in names(ometh)) {
    if (verbosity > 2) {
      cat("\t\t", meth, " ", sep="")
      pt <- proc.time()
    }
    res[[meth]] <-
      do.call(estimate, c(ometh[[meth]], opt_data, list(method = meth)))
    if (verbosity > 2) {
      cat((proc.time() - pt)[3], "\n")
    }
  }
  list(data=data, predict=res)
}

#' @export
create_nonparam_opt <- function(
  reps, n, sd, n_new,
  curve = c("spiral_closed", "spiral_open", "geodesic"), speed = pi,
  linfre = TRUE,
  lincos = TRUE,
  lingeo = TRUE,
  locfre = TRUE,
  locgeo = TRUE,
  trifre = TRUE,
  trigeo = 3) {

  o <- list(samp = list(),
            simu = list(),
            meth = list())

  o$simu$reps <- reps

  o$samp$n <- n
  o$samp$n_new <- n_new
  o$samp$noise <- "contracted"
  o$samp$sd <- sd

  curve <- match.arg(curve)
  switch(curve,
         spiral_closed = {
           o$samp$curve <- "spiral"
           o$samp$theta_min <- pi/4
           o$samp$theta_max <- pi/4
           o$samp$phi_start <- 0.5
           o$samp$circles <- 1
           periodic <- TRUE
         },
         spiral_open = {
           o$samp$curve <- "spiral"
           o$samp$theta_min <- pi/8
           o$samp$theta_max <- pi/8*7
           o$samp$phi_start <- 0.5
           o$samp$circles <- 1.5
           periodic <- FALSE
         },
         geodesic = {
           o$samp$curve <- "geodesic"
           o$samp$speed_bounds <- speed
         }
  )

  if (locfre) {
    o$meth$locfre$adapt <- "loocv"
    o$meth$locfre$kernel <- "epanechnikov"
    o$meth$locfre$bw <- 7
    o$meth$locfre$restarts <- 3
  }
  if (trifre) {
    o$meth$trifre$adapt <- "loocv"
    o$meth$trifre$num_basis <- 20
    o$meth$trifre$periodize <- !periodic
    o$meth$trifre$restarts <-  3
  }
  if (locgeo) {
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
  if (linfre) {
    o$meth$linfre$restarts <- 2
  }
  if (lincos) {
    o$meth$lincos$max_speed <- 10
    o$meth$lincos$restarts <- 2
  }
  if (lingeo) {
    o$meth$lingeo$max_speed <- 10
    o$meth$lingeo$restarts <- 2
  }
  o
}


#' @export
simulate <- function(opt_list, verbosity=1) {
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
      res[[j]] <- simu_run_one(opt$samp, opt$meth, verbosity)
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


#' Simulate using package parallel
#'
#' Repetitions of each setting a run in parallel.
#'
#' @export
simulate_parallel <- function(opt_list, cores=parallel::detectCores()) {
  cl <- parallel::makeCluster(cores)
  all_res <- lapply(opt_list, function(opt) {
    parallel::parLapply(
      cl,
      seq_len(opt$simu$reps),
      function(i) simu_run_one(osamp=opt$samp, ometh=opt$meth, verbosity=0))
  })
  parallel::stopCluster(cl)
  all_res
}

