# Data Generation for Uniform Example 

DGP <- function(mu = c(-6,0,6), n = 1000, N = 100, vmu = 2){
    histo <- function(y, x) 
	hist(y, x, plot = FALSE)$density
    x <- seq(-15, 15, len = 101)
    U <- matrix(runif(3 * N), N, 3)
    Y <- matrix(0, n, N)
    for(i in 1:N)
	Y[,i] <- sample(mu, n, prob = U[i,], replace = TRUE) + 
	    rnorm(1, sd = sqrt(2)) + rnorm(n)
    H <- apply(Y, 2, histo, x = x)
    H <- H/apply(H,2,sum)
    xm <- (x[-1] + x[-length(x)])/2
    M <- abs(outer(xm, xm, "-"))
    M <- M
    list(x = x, H = H, M = M)
}


require(Rmosek)
#set.seed(39)
n <- 100
Aq <- NA
bq <- NA
one <- matrix(1, 1, n)
A <- rbind(kronecker(one, Diagonal(n)),kronecker(Diagonal(n), one))
D <- DGP(N = 2)
x <- D$x
H <- D$H
f <- c(D$M)
b <- c(D$H)
lb <- rep(0, n^2)
ub <- rep(Inf, n^2)
P <- mosek_lptoprob(f,Aq,bq,A,b,lb,ub)
#r <- mosek(P, opts = list(verbose = 0))
r <- mosek(P) 
X <- matrix(r$sol$itr$xx,n,n)
xm <- (x[-1] + x[-length(x)])/2

X11(height = 6, width = 10)
par(mfrow = c(1,3))
for(i in 1:2){
    main <- paste("Histogram", i)
    plot(stepfun(x[-c(1,100)], H[,i]), main = main, do.points = FALSE)
}
G <- apply(H,2,cumsum)
xm <- (x[-1] + x[-length(x)])/2
plot(G[,1],xm, type = "l", xlab = "u", ylab = "Q(u)")
lines(G[,2],xm, col = 2)
g <- approxfun(xm, abs(G[,2] - G[,1]))
Ig <- integrate(g, -12, 12, subdivisions = 1000)

# What is the cost for random matching?
R <- 1000
X <- sample(xm, R, prob = H[,1], replace = TRUE)  
Y <- sample(xm, R, prob = H[,2], replace = TRUE)
MTC <- mean(abs(X - Y))

