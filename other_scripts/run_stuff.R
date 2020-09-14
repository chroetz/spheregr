
run_geodesic_reg <- function() {
  speed_max <- 10
  n <- 30
  data <- sample_data(n=n, sd=0.3, speed_max=speed_max)
  x <- data$x
  y <- data$y
  y_angle <- R32angle(y)
  x_true <- x
  y_true <- data$y_true
  y_true_angle <- R32angle(y_true)
  p_true <- data$p
  v_true <- data$v
  speed_true <- sqrt(sum(v_true^2))

  estim <- estimate_geodesic(x, y, x_true)
  estim_angle <- R32angle(estim)

  ise <- mean(dist_angle(y_true_angle, estim_angle)^2)
  cat("ISE:",ise, "\tspeed:", speed_true,"\n")

  plot(NA, pch=".",
       xlim=c(0, pi), ylim=c(0, 2*pi),
       xlab="alpha", ylab="phi", col="gray")
  striped_lines(y_true_angle, col=stripes_red, lwd=3)
  lines_jump(y_angle, col=rgb(0,0,0,0.3))
  points(y_angle)
  striped_lines(estim_angle, col=stripes_green, lwd=3)

  plot(NA, xlim=c(0, 1), ylim=c(0, 2*pi), xlab="t", ylab="phi")
  grid()
  lines_jump(cbind(x_true, y_true_angle[,2]), col=2, lwd=2)
  points(x, y_angle[,2])
  lines_jump(cbind(x_true, estim_angle[,2]), col=3, lwd=2)

  plot(NA, xlim=c(0, 1), ylim=c(0, pi), xlab="t", ylab="alpha")
  grid()
  lines_jump(cbind(x_true, y_true_angle[,1]), col=2, lwd=2)
  points(x, y_angle[,1])
  lines_jump(cbind(x_true, estim_angle[,1]), col=3, lwd=2)
}


run_frechet_reg <- function() {

  speed_max <- 10
  n <- 30
  data <- sample_data(n=n, sd=0.3, speed_max=speed_max)
  x <- data$x
  y <- data$y
  y_angle <- R32angle(y)
  x_true <- x
  y_true <- data$y_true
  y_true_angle <- R32angle(y_true)
  p_true <- data$p
  v_true <- data$v
  speed_true <- sqrt(sum(v_true^2))

  estim <- estimate_frechet(x, y, x_true)
  estim_angle <- R32angle(estim)

  ise <- mean(dist_angle(y_true_angle, estim_angle)^2)
  cat("ISE:",ise, "\tspeed:", speed_true,"\n")

  plot(NA, pch=".",
       xlim=c(0, pi), ylim=c(0, 2*pi),
       xlab="alpha", ylab="phi", col="gray")
  striped_lines(y_true_angle, col=stripes_red, lwd=3)
  lines_jump(y_angle, col=rgb(0,0,0,0.3))
  points(y_angle)
  striped_lines(estim_angle, col=stripes_green, lwd=3)

  plot(NA, xlim=c(0, 1), ylim=c(0, 2*pi), xlab="t", ylab="phi")
  grid()
  lines_jump(cbind(x_true, y_true_angle[,2]), col=2, lwd=2)
  points(x, y_angle[,2])
  lines_jump(cbind(x_true, estim_angle[,2]), col=3, lwd=2)

  plot(NA, xlim=c(0, 1), ylim=c(0, pi), xlab="t", ylab="alpha")
  grid()
  lines_jump(cbind(x_true, y_true_angle[,1]), col=2, lwd=2)
  points(x, y_angle[,1])
  lines_jump(cbind(x_true, estim_angle[,1]), col=3, lwd=2)
}

run_cosine_reg <- function() {
  speed_max <- 10
  n <- 30
  data <- sample_data(n=n, sd=0.3, speed_max=speed_max)
  x <- data$x
  y <- data$y
  y_angle <- R32angle(y)
  x_true <- x
  y_true <- data$y_true
  y_true_angle <- R32angle(y_true)
  p_true <- data$p
  v_true <- data$v
  speed_true <- sqrt(sum(v_true^2))

  estim <- estimate_cosine(x, y, x_true)
  estim_angle <- R32angle(estim)

  ise <- mean(dist_angle(y_true_angle, estim_angle)^2)
  cat("ISE:",ise, "\tspeed:", speed_true,"\n")

  plot(NA, pch=".",
       xlim=c(0, pi), ylim=c(0, 2*pi),
       xlab="alpha", ylab="phi", col="gray")
  striped_lines(y_true_angle, col=stripes_red, lwd=3)
  lines_jump(y_angle, col=rgb(0,0,0,0.3))
  points(y_angle)
  striped_lines(estim_angle, col=stripes_green, lwd=3)

  plot(NA, xlim=c(0, 1), ylim=c(0, 2*pi), xlab="t", ylab="phi")
  grid()
  lines_jump(cbind(x_true, y_true_angle[,2]), col=2, lwd=2)
  points(x, y_angle[,2])
  lines_jump(cbind(x_true, estim_angle[,2]), col=3, lwd=2)

  plot(NA, xlim=c(0, 1), ylim=c(0, pi), xlab="t", ylab="alpha")
  grid()
  lines_jump(cbind(x_true, y_true_angle[,1]), col=2, lwd=2)
  points(x, y_angle[,1])
  lines_jump(cbind(x_true, estim_angle[,1]), col=3, lwd=2)
}



run_linreg_reg <- function() {
  speed_max <- 10
  n <- 30
  data <- sample_data(n=n, sd=0.3, speed_max=speed_max)
  x <- data$x
  y <- data$y
  y_angle <- R32angle(y)
  x_true <- x
  y_true <- data$y_true
  y_true_angle <- R32angle(y_true)
  p_true <- data$p
  v_true <- data$v
  speed_true <- sqrt(sum(v_true^2))
  cat("speed:", speed_true,"\n")

  for (m in c("frechet", "geodesic", "cosine")) {
    cat("method:",m, "\t")
    estim <- linreg(x, y, x_true, method=m)
    estim_angle <- R32angle(estim)

    ise <- mean(dist_angle(y_true_angle, estim_angle)^2)
    cat("ISE:",ise, "\n")

    plot(NA, pch=".",
         xlim=c(0, pi), ylim=c(0, 2*pi),
         xlab="alpha", ylab="phi", col="gray")
    striped_lines(y_true_angle, col=stripes_red, lwd=3)
    lines_jump(y_angle, col=rgb(0,0,0,0.3))
    points(y_angle)
    striped_lines(estim_angle, col=stripes_green, lwd=3)

    plot(NA, xlim=c(0, 1), ylim=c(0, 2*pi), xlab="t", ylab="phi")
    grid()
    lines_jump(cbind(x_true, y_true_angle[,2]), col=2, lwd=2)
    points(x, y_angle[,2])
    lines_jump(cbind(x_true, estim_angle[,2]), col=3, lwd=2)

    plot(NA, xlim=c(0, 1), ylim=c(0, pi), xlab="t", ylab="alpha")
    grid()
    lines_jump(cbind(x_true, y_true_angle[,1]), col=2, lwd=2)
    points(x, y_angle[,1])
    lines_jump(cbind(x_true, estim_angle[,1]), col=3, lwd=2)
  }
}


