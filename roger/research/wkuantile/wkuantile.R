wquantile <- function(x, y, t = 0.5) {
#weighted quantile by brute force (full sorting) for scalar quantile regression
# through the origin model:  Q_{y_i} (t)  = x_i b
#
        ord <- order(y/x)
        b <- (y/x)[ord]
        wabs <- abs(x[ord])
        k <- sum(cumsum(wabs) < ((t - 0.5) * sum(x) + 0.5 * sum(wabs)))
        list(b = b[k + 1], k = ord[k + 1])
}

wkuantile <- function(x, y, tau = 0.5){
    n <- length(x)
    b <- y/x
    w <- abs(x)
    ip <- 1:n
    z <- .Fortran("wkuantile",
		  as.double(tau),
		  as.integer(n),
		  as.double(b),
		  as.double(w),
		  as.integer(ip),
		  q = double(1))
    z$q
}
dyn.load("wkuantile.so")
set.seed(17)
n = 1000001
y = rnorm(n)
x = runif(n)
tau = .8
t0 <- system.time(q0 <- wquantile(x,y,tau)$b)
t1 <- system.time(q1 <- wkuantile(x,y,tau))
