


run_linreg <- function(n=30, sd=0.3, speed_bounds = c(0, 10)) {
  p <- convert_a2e(c(2, 0))
  v <- matrix(c(0, 6, 0), nrow=1)
  data <- sample_data(n, sd, x_new=seq(0, 1, len=300), speed_bounds=speed_bounds,
                      p = p, v = v)
  cat("speed:", data$speed,"\n")

  res <- list()
  for (meth in linreg_methods) {
    cat("method:",meth, "\t")
    r <- linreg(data$x, data$y, data$x_new, method=meth, speed_bounds=speed_bounds)
    r$ise <- mean(dist(data$m, r$estim)^2)
    res[[meth]] <- r
  }

  plot(NA,
       xlim=c(0, pi), ylim=c(0, 2*pi),
       xlab="alpha", ylab="phi", col="gray")
  sphere_grid()
  striped_lines(data$m_a, col="red", lwd=3)
  #lines_jump(data$y_a, col=rgb(0,0,0,0.3))
  points(data$y_a)
  for (meth in all_methods) {
    striped_lines(res[[meth]]$estim_a, col=method_cols[[meth]], lwd=3)
    # points(res[[meth]]$estim_a, col=method_cols[[meth]])
  }

  plot(NA, xlim=c(0, 1), ylim=c(0, 2*pi), xlab="t", ylab="phi")
  grid()
  lines_jump(cbind(data$x_new, data$m_a[,2]), col=2, lwd=2)
  points(data$x, data$y_a[,2])
  for (meth in all_methods)
    lines_jump(cbind(data$x_new, res[[meth]]$estim_a[,2]), col=method_cols[[meth]], lwd=2)

  plot(NA, xlim=c(0, 1), ylim=c(0, pi), xlab="t", ylab="alpha")
  grid()
  lines_jump(cbind(data$x_new, data$m_a[,1]), col=2, lwd=2)
  points(data$x, data$y_a[,1])
  for (meth in all_methods)
    lines_jump(cbind(data$x_new, res[[meth]]$estim_a[,1]), col=method_cols[[meth]], lwd=2)
}



