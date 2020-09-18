load("sim_rand")

# mise
ise_stats <- sim_extract_ise_stats(sim_rand)
mise_table_from_ise_stats(ise_stats, c("mean", "sd"))

# speed
sim_extract_speed_mse(sim_rand)


load("sim_line")

# median curves

ise_stats <- sim_extract_ise_stats(sim_line)

idx <- 11
sim <- sim_line[[idx]]
idx_median <- ise_stats[[idx]]$idx_median
for (meth in linreg_methods)
  sim_plot_run(sim_line[[idx]][[idx_median[[meth]] ]])



# plot mean curves
plot(NA,
     xlim=c(0, 2*pi), ylim=c(0, pi),
     xlab="phi", ylab="alpha")
sphere_grid()
striped_lines(s[[1]]$m_a[,2:1], col="gray", lwd=3)
for (meth in linreg_methods) {
  mean_a <- mean_curves[[meth]][,1:2]
  geo_sd <- sqrt(mean_curves[[meth]][,3])
  striped_lines(mean_a[,2:1], col=method_colors[[meth]], lwd=3)
}
