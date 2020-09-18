striped_lines_a <- function(xy, ...) {
  striped_lines(xy[, 2:1], ...)
}

striped_lines <- function(xy, col, ...) {
  rgba <- col2rgb(col, alpha=TRUE)
  col_dark <- rgb(rgba[1]/2, rgba[2]/2, rgba[3]/2, rgba[4], maxColorValue=255)
  lines_jump(xy, col = col, lend = 1, ...)
  n <- nrow(xy)
  m <- 5
  for (i in 1:m) {
    lines_jump(
      xy[round(i * n / m - n / 2 / m):round(i * n / m), ],
      col = col_dark, lend = 1, ...)
  }
}

lines_jump <- function(y, ...) {
  diffs <- rowMeans(apply(y, 2, diff)^2)
  jumps <- which(diffs > 20 * median(diffs[diffs != 0]))
  segms <- c(0, jumps, nrow(y))
  for (i in 1:(length(segms) - 1)) {
    seg <- (segms[i] + 1):segms[i + 1]
    if (length(seg) < 2) next
    lines(y[seg, ], lend = 1, ...)
  }
}

sphere_grid <- function(n=7) {
  x <- matrix(seq(0, pi/2, len=100), ncol=1)
  for (phi in seq(0, 2*pi, len=n)) {
    p <- convert_a2e(c(pi/2, phi))

    v <- c(0, 0, 1)
    y <- Exp(p, x %*% v)
    y_a <- convert_e2a(y)
    lines(y_a[,2:1], col="gray")
    v <- c(0, 0, -1)
    y <- Exp(p, x %*% v)
    y_a <- convert_e2a(y)
    lines(y_a[,2:1], col="gray")
  }

  x <- matrix(seq(0, 2*pi, len=100), ncol=1)
  for (alpha in seq(pi/2/n, pi-pi/2/n, len=n)) {
    p <- convert_a2e(c(alpha, 0))

    v <- c(0, 1, 0)
    y <- Exp(p, x %*% v)
    y_a <- convert_e2a(y)
    lines(y_a[,2:1], col="gray")
  }
}


#' @export
sim_plot_run <- function(run) {
  plot(NA,
       xlim=c(0, 2*pi), ylim=c(0, pi),
       xlab="phi", ylab="alpha")
  sphere_grid()
  striped_lines_a(run$m_a, col="gray", lwd=3)
  points(run$y_a[, 2:1])
  for (meth in linreg_methods)
    striped_lines_a(run[[meth]]$estim_a, col=method_colors[[meth]], lwd=2)
  legend("topright", col=unlist(method_colors), lwd=2, legend=linreg_methods)
}


#' @export
sim_plot_biasvar <- function(sim) {
  estims_a <- array(NA, dim=c(
    samples = length(sim[[1]]$x_new),
    dims = 2,
    reps = length(sim))
  )
  mean_curves <- list()
  for (meth in linreg_methods) {
    for (i in seq_along(sim)) {
      estims_a[,,i] <- sim[[i]][[meth]]$estim_a
    }
    mean_curves[[meth]] <- t(apply(estims_a, 1, function(x) frechet_mean(t(x))))
  }

  x_new <- sim[[1]]$x_new
  m_a <- sim[[1]]$m_a

  par(mfrow=c(3,1), mar=c(0,4,1,1))
  plot(NA, xlim=c(0, 1), ylim=c(0, pi), xlab=NA, ylab="alpha", xaxt='n')
  grid()
  lines_jump(cbind(x_new, m_a[,1]), col="gray", lwd=3)
  for (meth in linreg_methods) {
    lines_jump(cbind(x_new, mean_curves[[meth]][,1]), col=method_colors[[meth]], lwd=3)
  }

  plot(NA, xlim=c(0, 1), ylim=c(0, 2*pi), xlab=NA, ylab="phi", xaxt='n')
  grid()
  lines_jump(cbind(x_new, m_a[,2]), col="gray", lwd=3)
  for (meth in linreg_methods) {
    lines_jump(cbind(x_new, mean_curves[[meth]][,2]), col=method_colors[[meth]], lwd=3)
  }

  par(mar=c(4,4,1,1))
  sd_max <- max(sapply(linreg_methods, function(meth) {
    max(sqrt(mean_curves[[meth]][,3]))
  }))
  plot(NA, xlim=c(0, 1), ylim=c(0, sd_max), xlab="x", ylab="sd")
  grid()
  for (meth in linreg_methods) {
    lines_jump(cbind(x_new, sqrt(mean_curves[[meth]][,3])), col=method_colors[[meth]], lwd=3)
  }
}

