# Figure 8.1 in van der Vaart -- Note that although this isn't mentioned in vdV 
# the MSEs in the figure are rescaled by sample size.  Note also that R's integrate
# function gets a little dicey for the n = 500 case, and I was too lazy to use the
# more reliable chi-squared route.

Risk <- function(t,n) {
    s <- sqrt(1/n)
    h <- function(x,t,s) x^2*dnorm(x, sd = s)
    R <- t
    for(i in 1:length(t)){ R[i] <- integrate(h,-Inf, -n^(-1/4) - t[i], s = s)$value +
        integrate(h,n^(-1/4) - t[i], Inf, s = s)$value +
        t[i]^2 * (pnorm(n^(-1/4)-t[i],sd = s) -  pnorm(-n^(-1/4)-t[i],sd = s) )
    }
    R
}
t <- -200:200/100
ns <-  c(500, 50, 5)
for(i in 1:length(ns)){
    n <- ns[i]
    if(i == 1) plot(t, n * Risk(t,n), col = i, type = "l")
    else lines(t, n * Risk(t,n), col = i)
}
