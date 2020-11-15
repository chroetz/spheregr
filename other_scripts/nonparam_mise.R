library(spheregr)
library(kableExtra)

reps <- 500
n_new <- 50
nonparam_methods <- c("locfre", "trifre", "locgeo", "trigeo")

opt_list <- list(
  create_opt(reps, 20, 0.25, n_new=n_new, curve="spiral_closed",  methods = nonparam_methods),
  create_opt(reps, 20, 0.25, n_new=n_new, curve="spiral_open",  methods = nonparam_methods),
  create_opt(reps, 80, 0.25, n_new=n_new, curve="spiral_closed",  methods = nonparam_methods),
  create_opt(reps, 80, 0.25, n_new=n_new, curve="spiral_open",  methods = nonparam_methods),
  create_opt(reps, 20, 1.00, n_new=n_new, curve="spiral_closed",  methods = nonparam_methods),
  create_opt(reps, 20, 1.00, n_new=n_new, curve="spiral_open",  methods = nonparam_methods),
  create_opt(reps, 80, 1.00, n_new=n_new, curve="spiral_closed",  methods = nonparam_methods),
  create_opt(reps, 80, 1.00, n_new=n_new, curve="spiral_open",  methods = nonparam_methods)
)

all_res <- simulate_parallel(opt_list, seed=1)

mise <- get_mise_tibble(opt_list, all_res)
pretty_print_mise(mise, "html")
