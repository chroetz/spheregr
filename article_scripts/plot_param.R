library(spheregr)

nn <- 1000
me <- c("linfre", "lingeo", "lincos")
cr <- "geodesic"
ac <- 0.5
gs <- 3

p <- matrix(c(1/sqrt(2),0,1/sqrt(2)), nrow=1)
v <- matrix(c(0,1,0), nrow=1)

opt_list <- list(
  create_opt(1, 10,  1/3, n_new=nn, curve=cr, geo_p=p, geo_v=v, geo_speed=3, methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 90,  1/3, n_new=nn, curve=cr, geo_p=p, geo_v=v, geo_speed=3, methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 10,  1.0, n_new=nn, curve=cr, geo_p=p, geo_v=v, geo_speed=3, methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 90,  1.0, n_new=nn, curve=cr, geo_p=p, geo_v=v, geo_speed=3, methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 10,  1/3, n_new=nn, curve=cr, geo_p=p, geo_v=v, geo_speed=6, methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 90,  1/3, n_new=nn, curve=cr, geo_p=p, geo_v=v, geo_speed=6, methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 10,  1.0, n_new=nn, curve=cr, geo_p=p, geo_v=v, geo_speed=6, methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 90,  1.0, n_new=nn, curve=cr, geo_p=p, geo_v=v, geo_speed=6, methods=me, accuracy=ac, grid_size=gs)
)

all_res <- simulate(opt_list, seed=1)

for (i in seq_along(all_res)) {
  file_name <- paste0("plot_param_", i, ".pdf")
  pdf(file_name, width=10, height=6)
  plot_run(all_res[[c(i,1)]], legend = "topright")
  title(main=paste0("n=", opt_list[[i]]$samp$n,", sd=",format(opt_list[[i]]$samp$sd, digits=2)),
        cex.main=1.5)
  dev.off()
}
