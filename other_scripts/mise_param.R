library(spheregr)

rp <- 1000
nn <- 100
me <- c("linfre", "lingeo", "lincos")
cr <- "geodesic"

opt_list <- list(
  create_opt(rp, 10,  0.1, n_new=nn, curve=cr, geo_speed=1.0, methods=me),
  create_opt(rp, 100, 0.1, n_new=nn, curve=cr, geo_speed=1.0, methods=me),
  create_opt(rp, 10,  1.0, n_new=nn, curve=cr, geo_speed=1.0, methods=me),
  create_opt(rp, 100, 1.0, n_new=nn, curve=cr, geo_speed=1.0, methods=me),
  create_opt(rp, 10,  0.1, n_new=nn, curve=cr, geo_speed=pi,  methods=me),
  create_opt(rp, 100, 0.1, n_new=nn, curve=cr, geo_speed=pi,  methods=me),
  create_opt(rp, 10,  1.0, n_new=nn, curve=cr, geo_speed=pi,  methods=me),
  create_opt(rp, 100, 1.0, n_new=nn, curve=cr, geo_speed=pi,  methods=me),
  create_opt(rp, 10,  0.1, n_new=nn, curve=cr, geo_speed=8.0, methods=me),
  create_opt(rp, 100, 0.1, n_new=nn, curve=cr, geo_speed=8.0, methods=me),
  create_opt(rp, 10,  1.0, n_new=nn, curve=cr, geo_speed=8.0, methods=me),
  create_opt(rp, 100, 1.0, n_new=nn, curve=cr, geo_speed=8.0, methods=me)
)

all_res <- simulate_parallel(opt_list, seed=1)

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
file_name <- paste0("simu_param_", timestamp, ".RData")
save(all_res, opt_list, file=file_name)

mise <- get_mise_tibble(opt_list, all_res, only_geo=TRUE)
pretty_print_mise(mise, "html")

