projected_mean <- function(y) {
  m <- norm_vec(colMeans(y))
}

var_ineq <- function(n=100, k=1e3, alpha=2, plot=FALSE, a=pi, b=2*pi) {
  y_a <- matrix(c(runif(n)*a, runif(n)*2*b), ncol=2)
  y <- convert_a2e(y_a)

  fm <- frechet_mean(y_a)
  m_in <- convert_a2e(fm[1:2])

  q_a <- matrix(c(runif(k)*pi, runif(k)*2*pi), ncol=2)
  q <- convert_a2e(q_a)
  y_dist <- apply(q, 1, function(qi) dist(matrix(rep(qi, each=nrow(y)), ncol=3), y))
  m_dist <- dist(q, m_in[rep(1, nrow(q)),])

  ex_d <- colMeans(y_dist^2)-fm[3]
  q_d <- m_dist^alpha

  if (plot) plot(q_d, ex_d/q_d, ylim = c(0,1), pch=".")

  min(ex_d/q_d)
}

res <- replicate(1000, var_ineq(alpha=1, k=1e3))
summary(res)




sd <- 1.0
n <- 1e5
k <- 1e3

y <- matrix(NA, nrow=n, ncol=3)
m <- matrix(c(1,0,0), ncol=3)
null <- pracma::nullspace(m)
for (i in 1:n) {
  noise_dir <- norm_vec(null %*% rnorm(2))
  y[i, ] <- Exp(m, t(noise_dir * rnorm(1, sd = sd)))
}
y_a <- convert_e2a(y)
fm <- frechet_mean(y_a)
m_in <- convert_a2e(fm[1:2])
sqrt(sum((m_in-m)^2))


q_a <- matrix(c(runif(k)*pi, runif(k)*2*pi), ncol=2)
q <- convert_a2e(q_a)
y_dist <- apply(q, 1, function(qi) dist(matrix(rep(qi, each=nrow(y)), ncol=3), y))
m_dist <- dist(q, m_in[rep(1, nrow(q)),])

ex_d <- colMeans(y_dist^2)-fm[3]
q_d <- m_dist^2

if (plot) plot(q_d, ex_d/q_d, ylim = c(0,1), pch=".")

min(ex_d/q_d)
