A <- function(x,r,s) {
  cos(x*s[1])*cos(x*s[2])*r[1] +
    cos(x*s[1])*sin(x*s[2])*r[2] +
    sin(x*s[1])*cos(x*s[2])*r[3] +
    sin(x*s[1])*sin(x*s[2])*r[4]
}
D <- function(x,r,s) {
  -sin(x*s[1])*sin(x*s[2])*r[1] +
    sin(x*s[1])*cos(x*s[2])*r[2] +
    cos(x*s[1])*sin(x*s[2])*r[3] -
    cos(x*s[1])*cos(x*s[2])*r[4]
}
B <- function(x,r,s) {
  (s[1]^2 + s[2]^2) * A(x,r,s) + 2*s[1]*s[2]*D(x,r,s)
}
C <- function(x,r,s) {
  s[1] *
    (sin(x*s[1])*cos(x*s[2])*r[1] +
       sin(x*s[1])*sin(x*s[2])*r[2] -
       cos(x*s[1])*cos(x*s[2])*r[3] -
       cos(x*s[1])*sin(x*s[2])*r[4]) +
    s[2] *
    (cos(x*s[1])*sin(x*s[2])*r[1] -
      cos(x*s[1])*cos(x*s[2])*r[2] +
      sin(x*s[1])*sin(x*s[2])*r[3] -
      sin(x*s[1])*cos(x*s[2])*r[4])
}
f0 <- function(x,r,s) {
  acos(A(x,r,s))
}
f2 <- function(x,r,s) {
  2 * acos(A(x,r,s)) * B(x,r,s) / (1 - A(x,r,s)^2)^(1/2) -
    2 * acos(A(x,r,s)) * A(x,r,s) * C(x,r,s)^2/(1-A(x,r,s)^2)^(3/2) +
    2* C(x,r,s)^2 / (1-A(x,r,s)^2)
}
r <- runif(4,min=-0.5,max=0.5)
s <- runif(2)*pi
x <- seq(-20,20,len=300)
y <- f2(x,r,s)
plot(x, y, type="l")
plot(x, f0(x,r,s), col=2, type="l")
grid()
max(A(x,r,s))
min(A(x,r,s))
