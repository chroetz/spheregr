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

sim_extract_ise_stats_one <- function(s) {
  opt <- list(n=s[[1]]$n, sd=s[[1]]$sd, speed=s[[1]]$speed)
  ise <- sim_extract_ise(s)
  summaries_ise <- lapply(ise, function(x) c(
    mean=mean(x), median=median(x), sd=sd(x), min=min(x), max=max(x),
    quartile1=quantile(x, 0.25), quartile3=quantile(x, 0.75)
  ))
  idx_median <- lapply(ise, which_median)
  median_curves <- list()
  for (meth in linreg_methods) {
    median_curves[[meth]] <- s[[idx_median[[meth]] ]][[meth]]$estim
  }
  list(
    opt = unlist(opt),
    summaries = summaries_ise,
    median_curves = median_curves,
    idx_median = idx_median
  )
}

#' @export
sim_extract_ise_stats <- function(sim) {
  lapply(sim, sim_extract_ise_stats_one)
}

#' @export
mise_table_from_ise_stats <- function(ise_stats, stats="mean") {
  extract_stat <- function(x, stat) {
    res <- sapply(linreg_methods, function(meth) x$summaries[[meth]][[stat]])
    names(res) <- paste0(linreg_methods, "_", stat, "_ise")
    res
  }
  extract_opts <- function(x) x$opt
  stat_tables <- lapply(
    stats,
    function(stat) sapply(ise_stats, extract_stat, stat=stat)
  )
  t(rbind(sapply(ise_stats, extract_opts), Reduce(rbind, stat_tables)))
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
  speed_sd <- lapply(speed_se, sd)
  c(opt, speed_mse=speed_mse, speed_sd=speed_sd)
}

#' @export
sim_extract_speed_mse <- function(sim) {
  lst <- lapply(sim, sim_extract_speed_mse_one)
  mat <- t(matrix(unlist(lst, recursive=TRUE), ncol=length(lst)))
  colnames(mat) <- c("n", "sd", "speed",
                     "geodesic_speed_mean_se", "cosine_speed_mean_se",
                     "geodesic_speed_sd_se", "cosine_speed_sd_se")
  tibble::as_tibble(mat)
}
