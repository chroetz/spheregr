striped_lines <- function(y, col, ...) {
  lines(y, col = col[1], lend = 1, ...)
  n <- nrow(y)
  m <- 5
  for (i in 1:(m - 0)) {
    lines(y[round(i * n / m - n / 2 / m):round(i * n / m), ], col = col[2], lend = 1, ...)
  }
}

lines_jump <- function(y, ...) {
  diffs <- rowMeans(apply(y, 2, diff)^2)
  jumps <- which(diffs > 100 * median(diffs[diffs != 0]))
  segms <- c(0, jumps, nrow(y))
  for (i in 1:(length(segms) - 1)) {
    lines(y[(segms[i] + 1):segms[i + 1], ], lend = 1, ...)
  }
}

stripes_red <- c(rgb(1, 0.5, 0.5), rgb(0.5, 0, 0))
stripes_green <- c(rgb(0.5, 1, 0.5), rgb(0, 0.5, 0))
stripes_blue <- c(rgb(0.5, 0.5, 1), rgb(0, 0, 0.5))
