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

