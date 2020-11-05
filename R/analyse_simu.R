get_ise <- function(all_res) {
  lapply(all_res, function(runs_res) {
    sapply(runs_res, function(res) {
      sapply(names(res$predict), function(meth) {
        mean(dist_a(res$predict[[meth]]$estim_a, res$data$m_new_a)^2)
      })
    })
  })
}

get_ise_summaries_from_ise <- function(ise) {
  lapply(ise, function(y) {apply(y, 1, function(x)
    c(
      mean=mean(x), median=median(x), sd=sd(x), min=min(x), max=max(x),
      quantile(x, 0.25), quantile(x, 0.75), idx_median=which_median(x)
    ))})
}

get_mise_tibble_from_ise_summaries <- function(opt_list, ise_summaries) {
  mise <- tibble::tibble()
  for (i in seq_along(opt_list)) {
    mise[i,"n"] <- opt_list[[i]]$samp$n
    mise[i,"sd"] <- opt_list[[i]]$samp$sd
    mise[i, "curve"] <- opt_list[[i]]$simu$curve
    mise[i, colnames(ise_summaries[[i]])] <- as.list(ise_summaries[[i]]["mean",])
  }
  mise
}

#' @export
get_mise_tibble <- function(opt_list, all_res) {
  ise <- get_ise(all_res)
  ise_summaries <- get_ise_summaries_from_ise(ise)
  get_mise_tibble_from_ise_summaries(opt_list, ise_summaries)
}


# TODO: all speed functions
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
