# target functions:
# * small circle
# *

small_circle <- function(n, theta_min, theta_max=theta_min, circles=1) {
  phi <- seq(0, 2*pi*circles, len=n) %% (2*pi)
  theta <- seq(theta_min, theta_max, len=n)
  m_a <- cbind(theta, phi)
  m_a
}


palette <- function(n) rainbow(n, s=1.0, v=0.8, alpha=NULL)
kernel <- function(x) 3/4*(1-x^2)*(abs(x)<1)
col_vec <- hsv((1:8 * (1 / ((1+sqrt(5))/2))) %% 1, s=1.0, v=0.9)
palette2 <- colorRampPalette(col_vec)
n <- 80
m_a <- small_circle(n, pi/8, pi/8*7, 1.5)
m <- convert_a2e(m_a)
y <- add_noise_contract(m, 0.5)
y_a <- convert_e2a(y)
h <- 0.2
x <- (1:n)/n
x_new <- seq(0, 1, len=20)

estim_a <- estimate_locgeo(x, y, x_new, kernel, h, speed_bounds=c(0,5))


plot_mercator_base()
lines_mercator(m_a, palette=palette, lwd=2)
points_mercator(y_a, palette=palette)
lines_mercator(estim_a, palette=palette, lwd=4)

