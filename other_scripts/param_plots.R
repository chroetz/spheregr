library(spheregr)

reps <- 1
n_new <- 1000
me <- c("linfre", "lingeo", "lincos")

p <- matrix(c(1/sqrt(2),0,1/sqrt(2)), nrow=1)
v <- matrix(c(0,1,0), nrow=1)

opt_list <- list(
  create_opt(reps, 10,  1.0, n_new=n_new, curve="geodesic", geo_p = p, geo_v = v, geo_speed = pi, methods = me),
  create_opt(reps, 20,  1.0, n_new=n_new, curve="geodesic", geo_p = p, geo_v = v, geo_speed = pi, methods = me),
  create_opt(reps, 40,  1.0, n_new=n_new, curve="geodesic", geo_p = p, geo_v = v, geo_speed = pi, methods = me),
  create_opt(reps, 80,  1.0, n_new=n_new, curve="geodesic", geo_p = p, geo_v = v, geo_speed = pi, methods = me),
  create_opt(reps, 100, 1.0, n_new=n_new, curve="geodesic", geo_p = p, geo_v = v, geo_speed = pi, methods = me)
)

all_res <- simulate(opt_list, seed=1)

for (res in all_res) {
  pdf(paste0("plots/plot_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".pdf"), width=10, height=6)
  plot_run(res[[1]])
  dev.off()
}
