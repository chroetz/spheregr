spiral <- function(n, theta_min, theta_max=theta_min, phi_start = 0.5, circles=1) {
  phi <- seq(phi_start, phi_start+2*pi*circles, len=n) %% (2*pi)
  theta <- seq(theta_min, theta_max, len=n)
  m_a <- cbind(theta, phi)
  m_a
}
palette <- function(n) rainbow(n, s=1.0, v=1.0, alpha=NULL)
palette_const <- function(v) function(n) rep(v, times=n)

ker_epan <- function(x) 3/4*(1-x^2)*(abs(x)<1)

n <- 60
contr <- 0.3

m_a <- spiral(n, pi/8, pi/8*7, circles=1.5)
m_new_a <- spiral(n_new, pi/8, pi/8*7, circles=1.5)
#m_a <- spiral(n, pi/4)
#m_new_a <- spiral(n_new, pi/4)

n_new <- 200
x_new <- seq(0, 1, len=n_new)
x <- seq(0, 1, len=n)
m <- convert_a2e(m_a)
#sd <- sqrt((pi^2-4)/2)*contr
y <- add_noise_contract(m, contr)
y_a <- convert_e2a(y)


cat("geo", system.time({
  geo <- estimate_geodesic(x, y, x_new, speed_bounds=c(0,5), restarts=7)
})[3], "\n")
geo_a <- geo$estim_a
geo_ise <- mean(dist_a(geo_a, m_new_a)^2)

cat("trigeo", system.time({
  trigeo <- estimate_trigeo(x, y, x_new, n_basis=5, speed_bounds=c(0,5))
})[3], "\n")#, "n_basis:", trigeo$n_basis,"\n")
trigeo_a <- trigeo$estim_a
trigeo_ise <- mean(dist_a(trigeo_a, m_new_a)^2)

cat("trifre", system.time({
  trifre <- estimate_trifre_loocv(x, y, x_new, periodize=TRUE)
})[3], "n_basis:", trifre$n_basis,"\n")
trifre_a <- trifre$estim_a
trifre_ise <- mean(dist_a(trifre_a, m_new_a)^2)

cat("locgeo", system.time({
  locgeo <- estimate_locgeo_loocv(x, y, x_new, kernel=ker_epan, speed_bounds=c(0,5))
})[3], "h:", locgeo$h, "\n")
locgeo_a <- locgeo$estim_a
locgeo_ise <- mean(dist_a(locgeo_a, m_new_a)^2)

cat("locfre", system.time({
  locfre <- estimate_locfre_loocv(x, y, x_new, kernel=ker_epan)
})[3], "h:", locfre$h, "\n")
locfre_a <- locfre$estim_a
locfre_ise <- mean(dist_a(locfre_a, m_new_a)^2)

print(data.frame(
  geo_ise,
  trifre_ise,
  locgeo_ise,
  locfre_ise))


pdf(paste0("plots/nonparam_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".pdf"), width=16, height=9)
layout(
  matrix(1:2, ncol=1),
  heights = c(9,1),
  widths = c(1)
)

plot_mercator_base()
lines_mercator(m_a, palette=palette_const(1), lwd=8)
lines_mercator(m_a, palette=palette_const(1), lwd=4)
lines_mercator(m_a, palette=palette, lwd=2)
points_mercator(y_a, m_a=m_a, palette=palette, pch=21, cex=1.5)

lines_mercator(geo_a, palette=palette_const(2), lwd=8)
lines_mercator(geo_a, palette=palette_const(1), lwd=4)
lines_mercator(geo_a, palette=palette, lwd=2)

lines_mercator(trifre_a, palette=palette_const(3), lwd=8)
lines_mercator(trifre_a, palette=palette_const(1), lwd=4)
lines_mercator(trifre_a, palette=palette, lwd=2)

lines_mercator(locgeo_a, palette=palette_const(4), lwd=8)
lines_mercator(locgeo_a, palette=palette_const(1), lwd=4)
lines_mercator(locgeo_a, palette=palette, lwd=2)

lines_mercator(locfre_a, palette=palette_const(5), lwd=8)
lines_mercator(locfre_a, palette=palette_const(1), lwd=4)
lines_mercator(locfre_a, palette=palette, lwd=2)

legend("bottomright", col=1:5, legend=c("true", "geo", "trifre", "locgeo", "locfre"), lwd=2)

par(mar=c(2,0,0,0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
x <- seq(0, 1, len=300)
colors=palette(300)
for (i in 1:299) rect(x[i], 0, x[i+1], 1, col=colors[i], border=NA)
ticks <- seq(0,1,0.2)
axis(1, at = ticks, labels = ticks, pos = 0)

dev.off()

