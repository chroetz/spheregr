
eval_trig_base <- function(x, num_basis) {
  k <- seq_len(num_basis/2)
  phi <- matrix(NA, nrow=2*length(k)+1, ncol=length(x))
  phi[1, ] <- 1
  phi[2*k, ] <- sqrt(2)*cospi(2 * k %*% t(x))
  phi[2*k+1, ] <- sqrt(2)*sinpi(2 * k %*% t(x))
  phi[1:num_basis, , drop=FALSE]
}


get_kernel_fun <- function(kernel = c("gaussian", "rectangular", "epanechnikov")) {
  if (is.character(kernel)) {
    kernel <- match.arg(kernel)
    kernel_fun <- switch(kernel,
                         gaussian = dnorm,
                         rectangular = function(x) as.numeric(abs(x) <= 0.5),
                         epanechnikov = function(x) 3 / 4 * (1 - x ^ 2) * (abs(x) < 1)
    )
  } else {
    kernel_fun <- kernel
  }
  kernel_fun
}


get_initial_parameters <- function(grid_size, max_speed=10, num_basis=1) {
  del <- 1/(2*grid_size)
  ticks <- seq(del, 1-del, len=grid_size)
  param_list <- c(
    list(alpha = pi * ticks, phi = 2*pi*ticks),
    rep(list(2*max_speed*ticks-max_speed), num_basis*2))
  as.matrix(expand.grid(param_list))
}
