library(spheregr)

rp <- 500
nn <- 50
me <- c("locfre", "trifre", "locgeo", "trigeo")

opt_list <- list(
  create_opt(rp, 20, 0.25, n_new=nn, curve="spiral_closed", methods=me),
  create_opt(rp, 20, 0.25, n_new=nn, curve="spiral_open",   methods=me),
  create_opt(rp, 80, 0.25, n_new=nn, curve="spiral_closed", methods=me),
  create_opt(rp, 80, 0.25, n_new=nn, curve="spiral_open",   methods=me),
  create_opt(rp, 20, 1.00, n_new=nn, curve="spiral_closed", methods=me),
  create_opt(rp, 20, 1.00, n_new=nn, curve="spiral_open",   methods=me),
  create_opt(rp, 80, 1.00, n_new=nn, curve="spiral_closed", methods=me),
  create_opt(rp, 80, 1.00, n_new=nn, curve="spiral_open",   methods=me)
)

all_res <- simulate_parallel(opt_list, seed=1)

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
file_name <- paste0("simu_nonparam_", timestamp, ".RData")
save(all_res, opt_list, file=file_name)

mise <- get_mise_tibble(opt_list, all_res)
pretty_print_mise(mise, "html")
