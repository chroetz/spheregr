#' @param p nx3
#' @param q nx3
#' @return nx3
Log <- function(q, p) {
  qTp <- rowSums(q * p) # nx1
  u <- p - qTp * q # nx3
  acos(crop1(qTp)) * u / row_norm3(u)
}

#' @param p 1x3
#' @param v nx3
#' @return nx3
Exp <- function(p, v) {
  nv <- row_norm3(v)
  cos(nv) %*% p + sin(nv) * v / nv
}

#' @param x,y nx2
#' @return vector length n
dist_angle <- function(x,y) {
  x <- matrix(x, ncol=2)
  y <- matrix(y, ncol=2)
  acos(cos(x[,1])*cos(y[,1])+sin(x[,1])*sin(y[,1])*cos(y[,2]-x[,2]))
}

#' @param x,y nx3
#' @return vector length n
dist <- function(x,y) {
  x <- matrix(x, ncol=3)
  y <- matrix(y, ncol=3)
  acos(rowSums(x*y))
}

#' @param x nx2
#' @return nx3
angle2R3 <- function(x) {
  x <- matrix(x, ncol=2)
  cbind(sin(x[,1])*cos(x[,2]),
        sin(x[,1])*sin(x[,2]),
        cos(x[,1]))
}

#' @param x nx3
#' @return nx2
R32angle <- function(x) {
  phi <- atan2(x[,2], x[,1])
  cbind(acos(x[,3]), ifelse(phi >= 0, phi, phi+2*pi))
}
