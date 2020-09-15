sim_extract_ise <- function(s) {
  ise_lst <- list(frechet=NA, geodesic=NA, cosine=NA)

  ise <- function(a, b) {
    mean(dist(a, b)^2)
  }

  for (i in seq_along(s)) {
    run <- s[[i]]
    for (meth in linreg_methods) {
      ise_lst[[meth]][i] <- ise(run$m, run[[meth]]$estim)
    }
  }

  ise_lst
}

sim_extract_mise_one <- function(s) {
  opt <- list(n=s[[1]]$n, sd=s[[1]]$sd, speed=s[[1]]$speed)
  ise <- sim_extract_ise(s)
  mise <- lapply(ise, mean)
  c(opt, mise)
}

sim_extract_mise <- function(sim) {
  lst <- lapply(sim, sim_extract_mise_one)
  mat <- t(matrix(unlist(lst, recursive=TRUE), ncol=length(lst)))
  colnames(mat) <- c("n", "sd", "speed", linreg_methods)
  tibble::as_tibble(mat)
}


sim_extract_speed_se <- function(s) {
  speed_se <- list(geodesic=NA, cosine=NA)
  for (i in seq_along(s)) {
    run <- s[[i]]
    speed_se$cosine[i] <- (run$cosine$speed - run$speed)^2
    speed_se$geodesic[i] <- (sqrt(sum(run$geodesic$v^2)) - run$speed)^2
  }
  speed_se
}

sim_extract_speed_mse_one <- function(s) {
  opt <- list(n=s[[1]]$n, sd=s[[1]]$sd, speed=s[[1]]$speed)
  speed_se <- sim_extract_speed_se(s)
  speed_mse <- lapply(speed_se, mean)
  c(opt, speed_se=speed_se)
}

sim_extract_speed_mse <- function(sim) {
  lst <- lapply(sim, sim_extract_speed_mse_one)
  mat <- t(matrix(unlist(lst, recursive=TRUE), ncol=length(lst)))
  colnames(mat) <- c("n", "sd", "speed", "cosine", "geodesic")
  tibble::as_tibble(mat)
}
