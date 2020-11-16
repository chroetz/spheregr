library(spheregr)

nn <- 1000
me <- c("locfre", "trifre", "locgeo", "trigeo")
ac <- 0.5
gs <- 3

opt_list <- list(
  create_opt(1, 10,  1/3, n_new=nn, curve="spiral_closed", methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 90,  1/3, n_new=nn, curve="spiral_closed", methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 10,  1.0, n_new=nn, curve="spiral_closed", methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 90,  1.0, n_new=nn, curve="spiral_closed", methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 10,  1/3, n_new=nn, curve="spiral_open",   methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 90,  1/3, n_new=nn, curve="spiral_open",   methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 10,  1.0, n_new=nn, curve="spiral_open",   methods=me, accuracy=ac, grid_size=gs),
  create_opt(1, 90,  1.0, n_new=nn, curve="spiral_open",   methods=me, accuracy=ac, grid_size=gs)
)

all_res <- simulate(opt_list, seed=1)

for (i in seq_along(all_res)) {
  file_name <- paste0("plot_nonparam_", i, ".pdf")
  pdf(file_name, width=10, height=6)
  plot_run(all_res[[c(i,1)]], legend = if (i %% 4 == 1) "topright")
  title(main=paste0("n=", opt_list[[i]]$samp$n,", sd=",format(opt_list[[i]]$samp$sd, digits=2)))
  dev.off()
}
