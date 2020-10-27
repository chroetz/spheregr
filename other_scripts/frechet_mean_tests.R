eps <- sqrt(.Machine$double.eps)/2
x <- seq(0+eps, 2*pi-eps, len=400)
b <- seq(-pi/2+eps, pi/2-eps, len=400)

energy <- function(a) {
  delta <- rowMeans((acos(sin(a)*sin(b) - cos(a)*cos(b)%*%t(cos(x))))^2)
  mean(delta * cos(pi/2-1.2*(pi/2+b)))
}

a <- seq(-pi/2, pi/2, len=100)
en <- sapply(a, energy)
plot(a, en, ylim=c(min(en), max(max(en), min(en)+0.5)), type="l")
grid()
which.min(en)


a <- seq(-pi/2, pi/2, len=400)
del <- function(a) {
  rowMeans((acos(sin(a)*sin(b) - cos(a)*cos(b)%*%t(cos(x))))^2)
}
d <- sapply(a, del)
plot(d %*% cos(a))
d2 <- apply(d*cos(a), 2, cumsum)
image(d2)

plot(d[1,], type="l", col=2)
lines(d[100,], type="l", col=3)
lines(d[150,], type="l", col=4)
lines(d[199,], type="l", col=5)
grid()

plot(d[,1], type="l", col=2)
lines(d[,100], type="l", col=3)
lines(d[,150], type="l", col=4)
lines(d[,199], type="l", col=5)
grid()

colMeans(d*cos(a))

d2 <- apply(d*cos(a), 2, cumsum)

image(d, asp=1, useRaster=TRUE)
image(d*cos(a), asp=1, useRaster=TRUE)
image(sqrt(d), asp=1, useRaster=TRUE)
plot(d[matrix(1:400,ncol=2, nrow=400)], type="l", col=2, lwd=2)
plot(d[matrix(c(1:400, 400:1),ncol=2, nrow=400)], type="l", col=2, lwd=2)
lines(2.65*a^2+3.32)



x <- seq(-1, 1, len=100)
y <- seq(-1, 1, len=100)
xy <- expand.grid(x=x,y=y)
f <- function(x, y) (acos(x)^2-acos(y)^2)/(y-x)
res_a <- f(xy[xy$y>xy$x,1], xy[xy$y>xy$x,2])
res_b <- f(xy[xy$y<xy$x,1], xy[xy$y<xy$x,2])
min(res_a)
max(res_b)
