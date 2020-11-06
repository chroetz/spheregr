library(spheregr)
set.seed(1)
reps <- 500
n_new <- 50
opt_list <- list(
  create_nonparam_opt(reps, 20, 0.25, n_new=n_new, curve="simple", trigeo=3),
  create_nonparam_opt(reps, 20, 0.25, n_new=n_new, curve="spiral", trigeo=3),
  create_nonparam_opt(reps, 80, 0.25, n_new=n_new, curve="simple", trigeo=3),
  create_nonparam_opt(reps, 80, 0.25, n_new=n_new, curve="spiral", trigeo=3),
  create_nonparam_opt(reps, 20, 1.00, n_new=n_new, curve="simple", trigeo=3),
  create_nonparam_opt(reps, 20, 1.00, n_new=n_new, curve="spiral", trigeo=3),
  create_nonparam_opt(reps, 80, 1.00, n_new=n_new, curve="simple", trigeo=3),
  create_nonparam_opt(reps, 80, 1.00, n_new=n_new, curve="spiral", trigeo=3)
) # reps=200, 45257.56
# reps=500: 103294.8

# cat("simulate", system.time({
#   all_res <- nonparam_simulate(opt_list, verbosity=10)
# })[3], "\n")
cat("simulate", system.time({
  all_res <- nonparam_simulate_parallel(opt_list)
})[3], "\n")

save(all_res, opt_list,
     file=paste0("simu/nonparam_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".RData"))

# pdf(paste0("plots/nonparam_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".pdf"), width=16, height=9)
# plot_run(all_res[[c(1,1)]])
# dev.off()

ise <- lapply(all_res, function(runs_res) {
  sapply(runs_res, function(res) {
    sapply(names(res$predict), function(meth) {
      mean(dist_a(res$predict[[meth]]$estim_a, res$data$m_new_a)^2)
    })
  })
})

mise <- lapply(ise, rowMeans)

library(tibble)
tbl <- tibble()
for (i in seq_along(opt_list)) {
  tbl[i,"n"] <- opt_list[[i]]$samp$n
  tbl[i,"sd"] <- opt_list[[i]]$samp$sd
  tbl[i, "curve"] <- opt_list[[i]]$simu$curve
  tbl[i, names(mise[[i]])] <- as.list(mise[[i]])
}

