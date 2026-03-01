# Simulation exercise for Q2 of 2012 574 exam

# Warning:  the following function shouldn't be evaluated at lambda == 0
# This is safe for the simulation below with even m.
lam <- function(lambda,x) (x^lambda - 1)/lambda
enorm <- function(x) sqrt(sum(x^2))
Gam <- function(lambda,x){
	Gam <- outer(lambda,x,lam)
	Gam/sqrt(apply(Gam^2,1,sum)) 
	}
gdotnorm <- function(x, m = 500){
	lambda <- seq(-1,1,length=m)
	G <- Gam(lambda,x)
	dG <- G[-1,] - G[-m,]
	sum(sqrt(apply(dG^2,1,sum)))
	}
phot <- function(w,x,m = 500, method = "JJ") {
# Two equivalent ways to get probability for "caps"
# The default seems to be the standard version as in my "Notes" 
# and the Johansen and Johnstone paper cited there.  The
# Student t version is from Knowles and Siegmund (ISR, 1989)
	n <- length(x)
	kappa <- gdotnorm(x, m = m) 
	if(method == "KS")
	return(kappa * ((1 - w^2)^((n-2)/2))/(2*pi) +
		(1 - pt(w * sqrt((n-1)/(1-w^2)), n-1)))
	else
	return(kappa * ((1 - w^2)^((n-2)/2))/(2*pi) +
		0.5 * (1 - pbeta(w^2, 0.5, (n-1)/2)))
	}
critval <- function(x, alpha = 0.05){
# Hotelling tube critical value for simple Box-Cox test 
   tube <- function(w, x, alpha) phot(w, x) -  alpha 
   uniroot(tube,c(0.05,0.5),x = x, alpha = alpha)$root
   }
Whot <- function(x,y,m = 500){
	y <- y/enorm(y)
	lambda <- seq(-1,1,length = m)
	G <- Gam(lambda,x)
	#plot(lambda, G %*%y, type = "l")
	max(G %*% y)
	}
BoxCox <- function(x,lambda) {
	if(lambda == 0) return(log(x)) 
	else (x^lambda - 1)/lambda
	}

require(doMC)
registerDoMC(8)
date()
sessionInfo()
set.seed(1939)
R <- 1000
lambdas <- c(-0.5,0,0.5)
cvs <- ns <- c(20, 50, 100, 500, 1000)
betas <- c(0,1,2)
xs <- list()
for(i in 1:length(ns)) xs[[i]] <- exp(rnorm(ns[i])) # one x vector per n
for(i in 1:length(ns)) cvs[i] <- critval(xs[[i]])
A <- array(0,c(length(ns),length(betas),R))
W <- rep(0,R)
ptime <- system.time({
   AA <- foreach(i = 1:length(lambdas)) %dopar% {
      for(j in 1:length(ns)){
         x <- xs[[j]]
         n <- ns[j]
         for(k in 1:length(betas)){
            for(r in 1:R){
               y <-  betas[k]/sqrt(n) * BoxCox(x,lambdas[i]) + rnorm(n)
               A[j,k,r] <- Whot(x,y)
               }
            }
         }
      A
    }
})
