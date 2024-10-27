## The huge speed up is because we have avoided all those multiplications
## by zero that contribute nothing anyway.
# nolint start
## 5.2
## 1
set.seed(5)
n <- 2000
w <- runif(n)
A <- matrix(runif(n * n), n, n)
system.time(B <- diag(w) %*% A)

## Let W = diag(w). Generically B_{ij} = \sum_k W_{ik} B_{kj}
## But W_{ik} = 0 except for W_{ii} = w_i, so B_{ij} = w_i B_{ij}
## given that B is actually stored column-wise that is in the order
## B_11 B_21, B_31, ... B_12 B_22 ... then the product can be
## computed using the recycling rule an vector multiplication...

system.time(B1 <- w * A)
range(B1 - B)
## The huge speed up is because we have avoided all those multiplications
## by zero that contribute nothing anyway.

## 2

julian(as.Date("2024-08-22"), origin = as.Date("2023-12-31"))

## 3
a <- factor(sample(c("fred", "george", "sue", "ann"), 20, replace = TRUE))
a
b <- factor(a, levels = c("ann", "sue", "fred", "george"))
b
as.numeric(a)
as.numeric(b)
## The codes are 1 for the first level, 2 for the second and so on
## so the change results in 4 -> 2, 2 -> 3 and 3 -> 4

al <- levels(a)
an <- as.numeric(a)
al[an]


## 5.4
x <- rnorm(20)
z <- rnorm(20)
y <- rnorm(20)
ii <- z != 0 & (z < 1 | y / z < 0)
x[ii] <- x[ii]^2

## 5.5 other options are possible...
## 1
x <- rnorm(20)
z <- rnorm(20)
y <- rnorm(20)
ii <- z != 0 & (z < 1 | y / z < 0)
x[ii] <- x[ii]^2

## 2
f <- function(x) log(x) - 1 / x

xr <- c(.1, 10) ## current bracketing interval
xt <- mean(xr) ## mid-point of interval
ft <- 1
while (abs(ft) > 1e-8) {
  xt <- mean(xr) ## mid-point of interval
  ft <- f(xt)
  if (ft > 0) xr[2] <- xt else xr[1] <- xt
}
xt
## 3
p <- 1000
system.time({
  A <- matrix(0, p, p)
  for (i in 1:p) for (j in i:p) A[i, j] <- 1
})

## 4
system.time(
  A <- matrix(rep(rep(1:0, p), rep(c(0, p), p) + rep(1:p, each = 2) * rep(c(1, -1), p)), p, p)
)


i <- rep(p, 2 * p)

i[seq(1, 2 * p, by = 2)] <- 1:p
i[seq(2, 2 * p, by = 2)] <- (p - 1):0
A <- matrix(rep(rep(1:0, p), i), p, p)

## 5
p <- 5
B <- matrix(1:(p * p), p, p)
A <- matrix(rep(rep(1:0, p), rep(c(0, p), p) + rep(1:p, each = 2) * rep(c(1, -1), p)), p, p)
B * A

## 5.6 mostly as above...

f <- function(x) log(x) - 1 / x
g <- function(x, a = 1, b = 1) x^b - exp(a * x) + 2

root <- function(xr, f, ...) {
  ## finds the solution of f(x,...) = 0 in the interval xr
  f0 <- f(xr[1], ...)
  f1 <- f(xr[2], ...)
  if (f0 == 0) {
    return(xr[1])
  }
  if (f1 == 0) {
    return(xr[2])
  }
  if (f0 * f1 > 0) stop("Interval does not appear to bracket solution")
  up <- sign(f1) ## +1 if function increasing, -1 otherwise
  repeat { ## bisection loop
    xt <- mean(xr) ## mid-point of interval
    ft <- f(xt, ...)
    if (up * ft > 0) xr[2] <- xt else xr[1] <- xt
    if (abs(ft) < 1e-8) break
  }
  xt
}

root(c(.1, 10), f)
root(c(0, 10), g, a = 1, b = 2)

# nolint end
