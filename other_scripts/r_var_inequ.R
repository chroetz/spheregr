
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

  if (plot) {
    plot(q_d, ex_d/q_d, ylim = c(0,1), pch=".")
    grid()
    abline(1,0, col=2)
    abline(0,0, col=2)
  }

 range(ex_d/q_d)
}

res <- replicate(100, var_ineq(n=5, k=1e2))
summary(res)
min(res[1,])
max(res[2,])
