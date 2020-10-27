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
  res <- cos(nv) %*% p + sin(nv) * v / nv
  b <- abs(nv) < 1e-8
  res[b,] <- p
  res
}

#' @param x,y nx2
#' @return vector length n
dist_a <- function(x, y) {
  x <- matrix(x, ncol = 2)
  y <- matrix(y, ncol = 2)
  acos(cos(x[, 1]) * cos(y[, 1]) + sin(x[, 1]) * sin(y[, 1]) * cos(y[, 2] - x[, 2]))
}

#' @param x,y nx3
#' @return vector length n
dist <- function(x, y) {
  x <- matrix(x, ncol = 3)
  y <- matrix(y, ncol = 3)
  acos(rowSums(x * y))
}

#' Convert angle representation to R3
#'
#' The first angle, theta, is in [0, pi] an describes longitude / the third
#' coordinate. The second angle, phi, is in [0, pi) an describes latitude /
#' rotation around the third axis.
#'
#' @param x nx2
#' @return nx3
convert_a2e <- function(x) {
  x <- matrix(x, ncol = 2)
  cbind(
    sin(x[, 1]) * cos(x[, 2]),
    sin(x[, 1]) * sin(x[, 2]),
    cos(x[, 1])
  )
}

#' @param x nx3
#' @return nx2
convert_e2a <- function(x) {
  phi <- atan2(x[, 2], x[, 1])
  cbind(acos(x[, 3]), ifelse(phi >= 0, phi, phi + 2 * pi))
}

#' Rotation matrix.
#'
#' A matrix that rotates (0,0,1) to m_a.
#' Note that there a many such matrices.
#' This matrix is the combination of a rotation around y- and z-axis.
#'
#' @param m_a Vector of length two. The two angles: theta and phi.
#' @return 3x3
rotation_matrix <- function(m_a) {
  matrix(c(
      cos(m_a[2])*cos(m_a[1]),
      sin(m_a[2])*cos(m_a[1]),
      -sin(m_a[1]),
      -sin(m_a[2]),
      cos(m_a[2]),
      0,
      cos(m_a[2])*sin(m_a[1]),
      sin(m_a[2])*sin(m_a[1]),
      cos(m_a[1])
    ),
    ncol=3,
    nrow=3
  )
}

#' Sample from the uniform distribution on the sphere.
#'
#' @param n number of samples
r_sphere_unif <- function(n) {
  u <- runif(n)
  theta <- acos(1-2*u)
  phi <- runif(n) * 2*pi
  cbind(theta, phi)
}

#' Sample from the contracted uniform distribution on the sphere.
#'
#' @param n number of samples
#' @param t double in [0,1]. contraction parameter. t=0 means constant m, t=1 is uniform distribution
#' @param m_a double vector of length 2. Point to contract to.
r_sphere_contract <- function(n, t, m_a=NULL, m=NULL) {
  x <- r_sphere_unif(n)
  x[,1] <- t*x[,1]
  if (is.null(m_a)) m_a <- convert_e2a(matrix(m, ncol=3))
  R <- rotation_matrix(m_a)
  convert_a2e(x) %*% t(R)
}

#' Add wrapped normal noise to points on sphere.
#'
#' @param m nx3. points on sohere in R3 coordinates.
#' @param sd positive double. standard deviation of normal distribution.
add_noise_normal <- function(m, sd) {
  for (i in seq_len(nrow(m))) {
    noise_dir <- norm_vec(t(pracma::nullspace(m[i, , drop = FALSE]) %*% rnorm(2)))
    m[i, ] <- Exp(m[i, ], noise_dir * rnorm(1, sd = sd))
  }
  m
}


#' Add contracted uniform noise to points on sphere.
#'
#' @param m nx3. points on sohere in R3 coordinates.
#' @param t double in [0,1]. Contraction parameter.
add_noise_contract <- function(m, t) {
  for (i in seq_len(nrow(m))) {
    m[i, ] <- r_sphere_contract(1, t, m=m[i, ])
  }
  m
}
