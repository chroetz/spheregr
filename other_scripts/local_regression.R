spiral <- function(n, theta_min, theta_max=theta_min, circles=1) {
  phi <- seq(0, 2*pi*circles, len=n) %% (2*pi)
  theta <- seq(theta_min, theta_max, len=n)
  m_a <- cbind(theta, phi)
  m_a
}
palette <- function(n) rainbow(n, s=1.0, v=1.0, alpha=NULL)
palette_const <- function(v) function(n) rep(v, times=n)

ker_epan <- function(x) 3/4*(1-x^2)*(abs(x)<1)

n <- 40
x <- seq(0, 1, len=n)
#m_a <- spiral(n, pi/8, pi/8*7, 1.5)
m_a <- spiral(n, pi/4)
#m_a <- cbind(
#  pi/2 + sin(seq(0, 2*pi, len=n)) * (pi/2-0.5),
#  seq(0, 2*pi, len=n))
m <- convert_a2e(m_a)
contraction <- 0.25
sd <- sqrt((pi^2-4)/2)*(1-contraction)
y <- add_noise_contract(m, 1-contraction)
y_a <- convert_e2a(y)

n_new <- 70
x_new <- seq(0, 1, len=n_new)
m_new_a <- spiral(n_new, pi/4)




plot_mercator_base()
lines_mercator(m_a, palette=palette_const(1), lwd=8)
lines_mercator(m_a, palette=palette_const(1), lwd=4)
lines_mercator(m_a, palette=palette, lwd=2)
points_mercator(y_a, m_a=m_a, palette=palette, pch=21)

geo <- estimate_geodesic(x, y, x_new, speed_bounds=c(0,5), restarts=7)$estim_a
lines_mercator(geo, palette=palette_const(2), lwd=8)
lines_mercator(geo, palette=palette_const(1), lwd=4)
lines_mercator(geo, palette=palette, lwd=2)

trifre <- estimate_trifre(x, y, x_new, n_basis=24)$estim_a
lines_mercator(trifre, palette=palette_const(3), lwd=8)
lines_mercator(trifre, palette=palette_const(1), lwd=4)
lines_mercator(trifre, palette=palette, lwd=2)

trifre <- estimate_trifre(x, y, x_new, n_basis=5)$estim_a
lines_mercator(trifre, palette=palette_const(4), lwd=8)
lines_mercator(trifre, palette=palette_const(1), lwd=4)
lines_mercator(trifre, palette=palette, lwd=2)

locgeo <- estimate_locgeo_loocv(x, y, x_new, ker_epan, speed_bounds=c(0,5))
locgeo_a <- locgeo$estim_a
ise <- mean(dist_a(locgeo_a, m_new_a)^2)
lines_mercator(locgeo_a, palette=palette_const(5), lwd=8)
lines_mercator(locgeo_a, palette=palette_const(1), lwd=4)
lines_mercator(locgeo_a, palette=palette, lwd=2)


locgeo_a <- list()
locfre_a <- list()
h <- c(0.05, 0.15, 0.5, 1.5)
for (i in 1:4) {
  locgeo_a[[i]] <- estimate_locgeo(x, y, x_new, ker_epan, h[i], speed_bounds=c(0,5))$estim_a
  locfre_a[[i]] <- estimate_locfre(x, y, x_new, ker_epan, h[i])$estim_a
}

par(mfrow=c(2,2))
for (i in 1:4) {
  plot_mercator_base()
  lines_mercator(m_a, palette=palette_const(1), lwd=8)
  lines_mercator(m_a, palette=palette_const(1), lwd=4)
  lines_mercator(m_a, palette=palette, lwd=2)
  points_mercator(y_a, m_a=m_a, palette=palette, pch=21)
  lines_mercator(locgeo_a[[i]], palette=palette_const(2), lwd=8)
  lines_mercator(locgeo_a[[i]], palette=palette_const(1), lwd=4)
  lines_mercator(locgeo_a[[i]], palette=palette, lwd=2)
  lines_mercator(locfre_a[[i]], palette=palette_const(3), lwd=8)
  lines_mercator(locfre_a[[i]], palette=palette_const(1), lwd=4)
  lines_mercator(locfre_a[[i]], palette=palette, lwd=2)
  title(main = paste0("h = ",h[i]))
}
