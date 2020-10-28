#' Calculate the Euclidean norm of an nx3 matrix rowwise.
#'
#' @param x nx3 matrix
#' @return vector of length n
row_norm3 <- function(x) sqrt(.Internal(rowSums(x^2, nrow(x), 3, FALSE)))

#' Return the normed vector.
#'
#' @param v numeric vector
#' @return the normed vector
norm_vec <- function(v) v / sqrt(sum(v^2))
