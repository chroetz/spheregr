#' @param p 1x3
#' @param v 1x3
#' @param x nx1
#' @param y nx3
geo_energy_w <- function(p, v, x, y, w) {
  z <- x %*% v # nx3
  n <- nrow(z)
  norm_z <- sqrt(.rowSums(z^2, n, 3L, FALSE))
  co <- tcrossprod(y, p) * cos(norm_z) # nx1
  si <- .rowSums(y * z, n, 3L, FALSE) / norm_z * sin(norm_z) # nx1
  s <- crop1(co + si)
  weighted.mean(acos(s)^2, w) / 2
}

optim_fn_w <- function(par, x, y, w) {
  p <- convert_a2e_1(par[1:2])
  v <- get_v_rcpp(p, par[3:4])
  dim(p) <- c(1L, 3L)
  dim(v) <- c(1L, 3L)
  geo_energy_w(p, v, x, y, w)
}

exec_optim <- function(x, y, w, speed_bounds, restarts = 2) {
  initial_parameters <- get_initial_parameters(speed_bounds, restarts)
  res_lst <- list()
  for (i in seq_len(nrow(initial_parameters))) {
    res_lst[[i]] <- stats::optim(
      initial_parameters[i, ], optim_fn_w,
      gr = NULL,
      x = x, y = y, w = w,
      method = "L-BFGS-B",
      lower = c(0, 0, -speed_bounds[2], -speed_bounds[2]),
      upper = c(pi, 2 * pi, speed_bounds[2], speed_bounds[2])
    )
  }
  values <- sapply(res_lst, function(x) x$value)
  idx <- which.min(values)
  res <- res_lst[[idx]]
  p <- convert_a2e_1(res$par[1:2])
  v <- get_v_rcpp(p, res$par[3:4])
  dim(p) <- c(1L, 3L)
  dim(v) <- c(1L, 3L)
  list(p = p, v = v)
}



plot_mercator_base <- function() {
  plot(NA, ylim=c(0, pi), xlim=c(0, 2*pi), ylab="theta", xlab="phi")
  sphere_grid()
}
lines_mercator <- function(y_a=NULL, y=NULL, palette=rainbow, ...) {
  if (is.null(y_a)) y_a <- convert_e2a(y)
  colors <- palette(nrow(y_a)-1)
  for (i in 1:(nrow(y_a)-1)) {
    lines(y_a[i:(i+1),2:1], col=colors[i], ...)
  }
}
points_mercator <- function(y_a=NULL, y=NULL, palette=rainbow, ...) {
  if (is.null(y_a)) y_a <- convert_e2a(y)
  points(y_a[,2:1], col=palette(nrow(y_a)), ...)
}


# target functions:
# * small circle
# *

small_circle <- function(theta, n) {
  phi <- seq(2*pi/n, 2*pi, len=n)
  m_a <- cbind(theta, phi)
  m_a
}

n <- 40
m_a <- small_circle(pi/4, n)
m <- convert_a2e(m_a)
y <- add_noise_contract(m, 0.5)
y_a <- convert_e2a(y)

palette <- function(n) rainbow(n, s=1.0, v=0.8, alpha=NULL)

kernel <- function(x) 3/4*(1-x^2)*(abs(x)<1)
x_plt <- seq(-2, 2, len=100)
plot(x_plt, kernel(x_plt), type="l")

weight <- function(x,t,h) {
  1/h*kernel((x-t)/h)
}

h <- 0.1
x <- (1:n)/n

estimate_at_t <- function(t) {
  w <- weight(x,t,h)
  w <- w / sum(w)
  res <- exec_optim(x, y, w, speed_bounds=c(0,5))
  estim <- Exp(res$p, t*res$v)
  estim_a <- convert_e2a(estim)
  estim_a
}

res_a <- t(sapply(seq(0,1,len=100), estimate_at_t))

plot_mercator_base()
lines_mercator(m_a, palette=palette, lwd=2)
points_mercator(y_a, palette=palette)
lines_mercator(res_a, palette=palette, lwd=4)
