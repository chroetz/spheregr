library(spheregr)

reps <- 1000
n_new <- 100
param_methods <- c("linfre", "lingeo", "lincos")

opt_list <- list(
  create_opt(reps, 10,  0.1, n_new=n_new, curve="geodesic", geo_speed = 1.0, methods = param_methods),
  create_opt(reps, 100, 0.1, n_new=n_new, curve="geodesic", geo_speed = 1.0, methods = param_methods),
  create_opt(reps, 10,  1.0, n_new=n_new, curve="geodesic", geo_speed = 1.0, methods = param_methods),
  create_opt(reps, 100, 1.0, n_new=n_new, curve="geodesic", geo_speed = 1.0, methods = param_methods),
  create_opt(reps, 10,  0.1, n_new=n_new, curve="geodesic", geo_speed = pi,  methods = param_methods),
  create_opt(reps, 100, 0.1, n_new=n_new, curve="geodesic", geo_speed = pi,  methods = param_methods),
  create_opt(reps, 10,  1.0, n_new=n_new, curve="geodesic", geo_speed = pi,  methods = param_methods),
  create_opt(reps, 100, 1.0, n_new=n_new, curve="geodesic", geo_speed = pi,  methods = param_methods),
  create_opt(reps, 10,  0.1, n_new=n_new, curve="geodesic", geo_speed = 8.0, methods = param_methods),
  create_opt(reps, 100, 0.1, n_new=n_new, curve="geodesic", geo_speed = 8.0, methods = param_methods),
  create_opt(reps, 10,  1.0, n_new=n_new, curve="geodesic", geo_speed = 8.0, methods = param_methods),
  create_opt(reps, 100, 1.0, n_new=n_new, curve="geodesic", geo_speed = 8.0, methods = param_methods)
)

all_res <- simulate_parallel(opt_list, seed=1)

save(all_res, opt_list, file=paste0("simu_param_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".RData"))

mise <- get_mise_tibble(opt_list, all_res)
pretty_print_mise(mise, "html")

