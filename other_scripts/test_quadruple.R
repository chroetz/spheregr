
n <- 1e4
q_a <- matrix(c(runif(4*n)*pi, runif(4*n)*2*pi), ncol=2)
q <- convert_a2e(q_a)

quadruple_const <- function(q_a) {
  lhs <- dist_a(q_a[1,], q_a[3,])^2 -
    dist_a(q_a[1,], q_a[4,])^2 -
    dist_a(q_a[2,], q_a[3,])^2 +
    dist_a(q_a[2,], q_a[4,])^2
  rhs <- dist_a(q_a[1,], q_a[2,]) * dist_a(q_a[3,], q_a[4,])
  abs(lhs) / rhs
}

quadruple_const_optim <- function(q_a) {
  q_a <- matrix(q_a, ncol=2)
  lhs <- dist_a(q_a[1,], q_a[3,])^2 -
    dist_a(q_a[1,], q_a[4,])^2 -
    dist_a(q_a[2,], q_a[3,])^2 +
    dist_a(q_a[2,], q_a[4,])^2
  rhs <- dist_a(q_a[1,], q_a[2,]) #* dist_a(q_a[3,], q_a[4,])
  - abs(lhs) / (rhs + 1e-14)
}

res <- sapply(1:n, function(i) quadruple_const(q_a[(4*i-3):(4*i),]))
max(res)
i_max <- which.max(res)
q_a_max <- q_a[(4*i_max-3):(4*i_max),]

lhs <- dist_a(q_a_max[1,], q_a_max[3,])^2 -
  dist_a(q_a_max[1,], q_a_max[4,])^2 -
  dist_a(q_a_max[2,], q_a_max[3,])^2 +
  dist_a(q_a_max[2,], q_a_max[4,])^2
rhs <- dist_a(q_a_max[1,], q_a_max[2,]) * dist_a(q_a_max[3,], q_a_max[4,])
abs(lhs) / rhs

res <- optim(q_a_max, quadruple_const_optim, method = "L-BFGS-B",
      lower = c(0,0,0,0,0,0,0,0),
      upper = c(rep(pi, 4), rep(2*pi, 4)))
q_a_max <- res$par
lhs <- dist_a(q_a_max[1,], q_a_max[3,])^2 -
  dist_a(q_a_max[1,], q_a_max[4,])^2 -
  dist_a(q_a_max[2,], q_a_max[3,])^2 +
  dist_a(q_a_max[2,], q_a_max[4,])^2
rhs <- dist_a(q_a_max[1,], q_a_max[2,])
abs(lhs) / rhs
