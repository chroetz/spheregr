#' Calculate the Euclidean norm of an nx3 matrix rowwise.
#'
#' @param x nx3 matrix
#' @return vector of length n
row_norm3 <- function(x) sqrt(.Internal(rowSums(x^2, nrow(x), 3, FALSE)))

#' Set values larger than 1 to 1 and smaller than -1 to -1.
#'
#' @param x numeric vector
#' @return the vector cropped to [-1,1]
crop1 <- function(x) {
  x[x < -1] <- -1
  x[x > 1] <- 1
  x
}

#' Return the normed vector.
#'
#' @param v numeric vector
#' @return the normed vector
norm_vec <- function(v) v / sqrt(sum(v^2))
