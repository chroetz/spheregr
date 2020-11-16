which_median <- function(x) {
  n <- length(x)
  if (n %% 2 == 0) order(x)[n/2] else order(x)[n/2+1]
}
