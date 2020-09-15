striped_lines <- function(xy, col, ...) {
  col_dark <- rgb(t(col2rgb(col)*0.5), maxColorValue=255)
  lines_jump(xy, col = col, lend = 1, ...)
  n <- nrow(xy)
  m <- 5
  for (i in 1:(m - 0)) {
    lines_jump(xy[round(i * n / m - n / 2 / m):round(i * n / m), ], col = col_dark, lend = 1, ...)
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
