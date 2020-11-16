n <- 1000
colors <- rainbow(n, alpha=1)
pdf("plot_rainbow_scale.pdf", width=20, height=2)
plot(NA, xlim=c(0,1), ylim=c(0,1), yaxt='n', xlab="t", ylab=NA, bty="n")
for (i in 1:n)
  rect((i-1)/n, 0, i/n, 1, col=colors[i], border=NA)
dev.off()
