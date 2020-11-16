#' @param p 1x3
#' @param v kx3
#' @param phi nxk
#' @param y nx3
trigeo_energy <- function(p, v, phi, y) {
  z <- phi %*% v # nx3
  n <- nrow(z)
  norm_z <- sqrt(.rowSums(z^2, n, 3L, FALSE))
  co <- tcrossprod(y, p) * cos(norm_z) # nx1
  si <- .rowSums(y * z, n, 3L, FALSE) / norm_z * sin(norm_z) # nx1
  s <- clamp1(co + si)
  sum(acos(s)^2) / 2
}

trigeo_optim_fn <- function(par, phi, y, num_basis) {
  p <- convert_a2e_1(par[1:2])
  v <- matrix(NA, nrow=num_basis, ncol=3)
  for (j in 1:num_basis) v[j, ] <- get_v_rcpp(p, par[c(2*j+1, 2*j+2)])
  dim(p) <- c(1L, 3L)
  trigeo_energy(p, v, phi, y)
}

#' @param accuracy double in [0,1]. trade optimization accuracy for
#'   computational speed
trigeo_optim <- function(x, y, max_speed, num_basis, grid_size, accuracy) {
  init_par <- get_initial_parameters(grid_size, max_speed, num_basis=num_basis)
  res_lst <- list()
  phi <- t(eval_trig_base(x, num_basis+1))[,-1,drop=FALSE]
  for (i in seq_len(nrow(init_par))) {
    res_lst[[i]] <- stats::optim(
      init_par[i, ], trigeo_optim_fn,
      gr = NULL,
      phi = phi, y = y, num_basis = num_basis,
      method = "L-BFGS-B",
      lower = c(0, 0, rep(-max_speed, times=2*num_basis)),
      upper = c(pi, 2 * pi, rep(max_speed, times=2*num_basis)),
      control = list(factr = 10^(17 * (1 - accuracy)))
    )
  }
  values <- sapply(res_lst, function(x) x$value)
  idx <- which.min(values)
  res <- res_lst[[idx]]
  p <- convert_a2e_1(res$par[1:2])
  v <- matrix(NA, nrow=num_basis, ncol=3)
  for (j in 1:num_basis) v[j, ] <- get_v_rcpp(p, res$par[c(2*j+1, 2*j+2)])
  dim(p) <- c(1L, 3L)
  list(p = p, v = v, par=res$par)
}

.estimate_trigeo <- function(x, y, x_new, num_basis, periodize, max_speed, grid_size, accuracy) {

  if (periodize) {
    n <- length(x)
    x <- c(x, rev(2-x))/2
    y <- rbind(y, y[n:1,])
    x_new <- x_new/2
  }

  res <- trigeo_optim(x, y, max_speed, num_basis, grid_size, accuracy)
  phi <- t(eval_trig_base(x_new, num_basis+1))[,-1,drop=FALSE]
  estim <- Exp(res$p, phi %*% res$v)
  estim_a <- convert_e2a(estim)
  c(list(estim=estim, estim_a=estim_a), res)
}


#' @export
estimate_trigeo <- function(
  x, y, x_new,
  adapt=c("loocv", "none"), num_basis=5, periodize=FALSE,
  max_speed = 10, grid_size = 2, accuracy=0.3
) {
  adapt <- match.arg(adapt)
  if (adapt == "loocv") {
    nums_basis <- 1:num_basis
    dists <- sapply(nums_basis, function(num_b) {
      v <- sapply(seq_along(x), function(j) {
        res <- .estimate_trigeo(x[-j], y[-j,], x[j],
                                num_b, periodize, max_speed,
                                grid_size, accuracy)
        dist(res$estim, y[j, ])
      })
      mean(v)
    })
    num_b <- nums_basis[which.min(dists)]
  } else {
    num_b <- num_basis
  }
  c(.estimate_trigeo(x, y, x_new, num_b, periodize, max_speed, grid_size, accuracy),
    list(num_basis = num_b))
}
