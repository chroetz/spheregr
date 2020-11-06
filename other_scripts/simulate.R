library(spheregr)
set.seed(1)
reps <- 4
n_new <- 30
param_methods <- c("linfre", "lingeo", "lincos")
nonparam_methods <- c("locfre", "trifre", "locgeo", "trigeo")
all_methods <- c("linfre", "lincos", "lingeo", "locfre", "trifre", "locgeo", "trigeo")
opt_list <- list(
  create_opt(reps, 20, 0.25, n_new=n_new, curve="geodesic",  methods = param_methods),
  create_opt(reps, 40, 0.25, n_new=n_new, curve="geodesic",  methods = param_methods)
)

#all_res <- simulate(opt_list, verbosity=5)
all_res <- simulate_parallel(opt_list, seed=1)

mise <- get_mise_tibble(opt_list, all_res)

plot_run(all_res[[1]][[2]])
