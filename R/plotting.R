palette_rainbow <- function(n) rainbow(n, s=1.0, v=1.0, alpha=NULL)
palette_const <- function(v) function(n) rep(v, times=n)

method_colors <- c(
  linfre = "#FF7F00",
  lincos = "#00FFFF",
  lingeo = "#7F3F7F",
  locfre = "#FF0000",
  trifre = "#00FF00",
  locgeo = "#0000FF",
  trigeo = "#FF00FF"
)

method_uppercamelcase <- c(
  linfre = "LinFre",
  lincos = "LinCos",
  lingeo = "LInGeo",
  locfre = "LocFre",
  trifre = "TriFre",
  locgeo = "LocGeo",
  trigeo = "TriGeo"
)

plot_mercator_base <- function(axis_cex=1.0, latex=FALSE) {
  plot(NA, ylim=c(0, pi), xlim=c(0, 2*pi), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
  xtick <- seq(0, 2*pi, len=5)
  ytick <- seq(0, pi, len=3)
  if (latex) {
    xticktxt <- c("$0$", "$\\frac12\\pi$", "$\\pi$", "$\\frac32\\pi$", "$2\\pi$")
    yticktxt <- c("$0$", "$\\frac12\\pi$", "$\\pi$")
    title(xlab="$\\varphi$", ylab="$\\vartheta$")
  } else {
    xticktxt <- c("0", "pi/2", "pi", "3pi/2", "2pi")
    yticktxt <- c("0", "pi/2", "pi")
    title(xlab="phi", ylab="theta", cex.lab=axis_cex)
  }
  axis(side=1, at=xtick, labels = xticktxt, cex.axis=axis_cex)
  axis(side=2, at=ytick, labels = yticktxt, cex.axis=axis_cex)
  sphere_grid()
}

wrap_line <- function(theta, phi, ...) {
  if (phi[1] < pi/2 & phi[2] > 2*pi-pi/2) {
    lines(c(phi[1],phi[2]-2*pi), theta, ...)
    lines(c(phi[1]+2*pi,phi[2]), theta, ...)
  } else if (phi[2] < pi/2 & phi[1] > 2*pi-pi/2) {
    lines(c(phi[1],phi[2]+2*pi), theta, ...)
    lines(c(phi[1]-2*pi,phi[2]), theta, ...)
  } else {
    lines(phi, theta, ...)
  }
}

lines_mercator <- function(y_a=NULL, y=NULL, palette=rainbow, ...) {
  if (is.null(y_a)) y_a <- convert_e2a(y)
  colors <- palette(nrow(y_a)-1)
  for (i in 1:(nrow(y_a)-1)) {
    wrap_line(y_a[i:(i+1),1], y_a[i:(i+1),2], col=colors[i], ...)
  }
}

points_mercator <- function(y_a=NULL, y=NULL, m_a=NULL, palette=rainbow, ...) {
  if (is.null(y_a)) y_a <- convert_e2a(y)
  colors <- palette(nrow(y_a))
  points(y_a[,2:1], bg=colors, col="black", ...)
  if (!is.null(m_a)) {
    for (i in 1:nrow(m_a)) {
      wrap_line(c(y_a[i,1], m_a[i,1]), c(y_a[i,2], m_a[i,2]), col=colors[i], lwd=0.5)
    }
  }
}

#' @export
plot_run <- function(res, rainbow=TRUE, legend=FALSE, legend_lwd=8, legend_cex=1.5, axis_cex=1.5) {
  with(res, {
    with(data, {
      plot_mercator_base(axis_cex)
      if (rainbow) {
        lines_mercator(m_new_a, palette=palette_const(1), lwd=8)
        lines_mercator(m_new_a, palette=palette_rainbow, lwd=2)
        points_mercator(y_a, m_a=m_a, palette=palette_rainbow, pch=21, cex=1.5)
      } else {
        lines_mercator(m_new_a, palette=palette_const(1), lwd=2)
        points_mercator(y_a, m_a=m_a, palette=palette_const(1), pch=21, cex=1.5)
      }
    })
    for (meth in names(predict)) {
      with(predict[[meth]], {
        if (rainbow) {
          lines_mercator(estim_a, palette=palette_const(method_colors[meth]), lwd=8)
          lines_mercator(estim_a, palette=palette_const(1), lwd=4)
          lines_mercator(estim_a, palette=palette_rainbow, lwd=2)
        } else {
          lines_mercator(estim_a, palette=palette_const(method_colors[meth]), lwd=2)
        }
      })
    }
  })

  if (!is.null(legend) && !isFALSE(legend) && !is.na(legend)) {
    par_family <- par("family")
    par(family="mono")
    pos <- "topright"
    if (is.character(legend)) pos <- legend
    methods <- names(res$predict)
    legend(
      pos,
      col=c("black", method_colors[methods]),
      legend=c("true", method_uppercamelcase[methods]), lwd=legend_lwd, cex=legend_cex, bg="white")
    par(family = par_family)
  }
}

plot_rainbow_time <- function() {
  #mar <- par("mar")
  #par(mar=c(2,0,0,0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  x <- seq(0, 1, len=300)
  colors=palette_rainbow(300)
  for (i in 1:299) rect(x[i], 0, x[i+1], 1, col=colors[i], border=NA)
  ticks <- seq(0,1,0.2)
  axis(1, at = ticks, labels = ticks, pos = 0)
  #par(mar=mar)
}


lines_jump <- function(y, ...) {
  diffs <- rowMeans(apply(y, 2, diff)^2)
  jumps <- which(diffs > 20 * median(diffs[diffs != 0]))
  segms <- c(0, jumps, nrow(y))
  for (i in 1:(length(segms) - 1)) {
    seg <- (segms[i] + 1):segms[i + 1]
    if (length(seg) < 2) next
    lines(y[seg, ], lend = 1, ...)
  }
}

sphere_grid <- function(n=7) {
  x <- matrix(seq(0, pi/2, len=100), ncol=1)
  for (phi in seq(0, 2*pi, len=n)) {
    p <- convert_a2e(c(pi/2, phi))

    v <- c(0, 0, 1)
    y <- Exp(p, x %*% v)
    y_a <- convert_e2a(y)
    lines(y_a[,2:1], col="gray")
    v <- c(0, 0, -1)
    y <- Exp(p, x %*% v)
    y_a <- convert_e2a(y)
    lines(y_a[,2:1], col="gray")
  }

  x <- matrix(seq(0, 2*pi, len=100), ncol=1)
  for (alpha in seq(pi/2/n, pi-pi/2/n, len=n)) {
    p <- convert_a2e(c(alpha, 0))

    v <- c(0, 1, 0)
    y <- Exp(p, x %*% v)
    y_a <- convert_e2a(y)
    lines(y_a[,2:1], col="gray")
  }
}


# TODO

#' @export
sim_plot_biasvar <- function(sim) {
  estims_a <- array(NA, dim=c(
    samples = length(sim[[1]]$x_new),
    dims = 2,
    reps = length(sim))
  )
  mean_curves <- list()
  for (meth in linreg_methods) {
    for (i in seq_along(sim)) {
      estims_a[,,i] <- sim[[i]][[meth]]$estim_a
    }
    mean_curves[[meth]] <- t(apply(estims_a, 1, function(x) frechet_mean(t(x))))
  }

  x_new <- sim[[1]]$x_new
  m_a <- sim[[1]]$m_a

  par(mfrow=c(3,1), mar=c(0,4,1,1))
  plot(NA, xlim=c(0, 1), ylim=c(0, pi), xlab=NA, ylab="alpha", xaxt='n')
  grid()
  lines_jump(cbind(x_new, m_a[,1]), col="gray", lwd=3)
  for (meth in linreg_methods) {
    lines_jump(cbind(x_new, mean_curves[[meth]][,1]), col=method_colors[[meth]], lwd=3)
  }

  plot(NA, xlim=c(0, 1), ylim=c(0, 2*pi), xlab=NA, ylab="phi", xaxt='n')
  grid()
  lines_jump(cbind(x_new, m_a[,2]), col="gray", lwd=3)
  for (meth in linreg_methods) {
    lines_jump(cbind(x_new, mean_curves[[meth]][,2]), col=method_colors[[meth]], lwd=3)
  }

  par(mar=c(4,4,1,1))
  sd_max <- max(sapply(linreg_methods, function(meth) {
    max(sqrt(mean_curves[[meth]][,3]))
  }))
  plot(NA, xlim=c(0, 1), ylim=c(0, sd_max), xlab="x", ylab="sd")
  grid()
  for (meth in linreg_methods) {
    lines_jump(cbind(x_new, sqrt(mean_curves[[meth]][,3])), col=method_colors[[meth]], lwd=3)
  }
}

