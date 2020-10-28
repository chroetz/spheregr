crop1 <- function(x) {
  x[x < -1] <- -1
  x[x > 1] <- 1
  x
}

crop2 <- function(x) {
  ifelse(x > 1, 1, ifelse(x < -1, -1, x))
}

crop3 <- function(x) {
  pmin(1, pmax(-1, x))
}

Rcpp::cppFunction('
NumericVector crop_mm(NumericVector x) {
  for(int i = 0; i < x.size(); ++i) {
    x[i] = std::max(-1.0, std::min(x[i], 1.0));
  }
  return x;
}')

Rcpp::cppFunction('
NumericVector crop_if(NumericVector x) {
  for(int i = 0; i < x.size(); ++i) {
    if (x[i] < -1) x[i] = -1;
    if (x[i] > 1) x[i] = 1;
  }
  return x;
}')


n <- 1e6
x <- runif(n, min=-2, max=2)
bm <- bench::mark(
  crop_mm = crop_mm(x),
  crop_if = crop_if(x),
  iterations=1e3
)
print(bm)
plot(bm)
