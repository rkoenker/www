# Simulation a la Dunson and Taylor (2005)

# Jefferies "substitute likelihood" a la Lavine (1995)
LJeff <- function(x,p,Q,Big = .Machine$double.xmax){
	if(any(diff(Q) < 0)) return(0)
	n <- length(x)
	u <- hist(x,c(-Big,Q,Big),plot = FALSE)$counts
	dp <- diff(c(0,p,1))
	exp(lgamma(n+1) - sum(lgamma(u+1)) + sum(u * log(dp)))
	}
	
# Prior distribution Gaussian(mu,Sigma)

dmvnorm <- function(x, mu, SigmaInv, detFactor = NULL, K = 1/sqrt(2*pi)){
	p <- length(x)
	if(!length(detFactor)) detFactor = sqrt(det(SigmaInv))
	(K^p) * detFactor * exp( -.5 * t(x - mu) %*% SigmaInv %*% (x-mu))
	}
phi0 <- function(x, mu = rep(0,length(x)), SigmaInv = diag(length(x)),detFactor)
	dmvnorm(x,mu,SigmaInv,detFactor)


# Sample quantiles and their (asymptotic) covariance matrix
Qtile <- function(x,taus) {
   n <- length(x)
   fit <- lapply(summary(rq(x ~ 1, tau = taus),se = "iid"),coefficients)
   fit <- matrix(unlist(fit),4)
   b <- fit[1,]
   Vb <- fit[2,]^2
   Omega <- outer(taus,taus,"pmin") - outer(taus,taus) # Brownian Bridge 
   f <- sqrt(diag(Omega)/(n * Vb)) # Density estimates
   V <- Omega/(n * outer(f,f)) # Avar(b)
   list(Qhat = b, V = V)
   }

# Metropolis Loop
	
MetroCycle <- function(a, x, p, S, a0, S0inv, detFactor, kappa = 1){
	at <- mvrnorm(1, a, kappa * S)
	if(any(diff(at) < 0)) 
		at <- a
	else {
		r <- min(1,(LJeff(x,p,at) * phi0(at,detFactor = detFactor)) /
			(LJeff(x,p,a) * phi0(a,detFactor = detFactor)))
                if(is.nan(r)) at <- a
		else if(runif(1) > r) at <- a
		}
	at
	}

QMCMC <- function(x, taus, a0, S0, J = 5000, Burnin = 1000){
	# Initialization of Proposal Distribution
	S0inv <- solve(S0)
	detFactor <- 1/sqrt(det(S0))
	z <- Qtile(x,taus)
	a1 <- z$Qhat
	S1 <- z$V
	S1inv <- solve(S1)
	S <- solve(S0inv + S1inv)
	a <- S %*% (S0inv %*% a0 + S1inv %*% a1)

	# MCMC Loop
	A <- matrix(0,length(a),J)
	for(j in 1:J){
		a <- MetroCycle(a, x, taus, S, a0, S0inv, 
			detFactor = detFactor, kappa = 1)
		A[,j] <- a
		}
	A <- A[,-(1:Burnin)]
	Qhat <- apply(A,1,mean)
	Cred <- apply(A,1,quantile,probs = c(.025,.05,.95,.975))
	AR <- length(rle(A[1,])$lengths)/(J - Burnin)
	list(Qhat = Qhat, Cred = Cred, AcceptanceRatio = AR)
	}

require(quantreg) #for rq
require(MASS) # for mvrnorm
require(nor1mix) # for normal mixtures

# Three target densities as in Dunson and Taylor
nM2 <- norMix(mu = c(-0.6,-0.5,-0.4,0.5,1),sig2 = c(0.4,0.5,0.6,0.95,1.15)^2)
nM3 <- norMix(mu = c(-0.5,-0.5,-0.5,0.75,1.5),sig2 = c(0.25,0.25,0.25,0.5,1)^2)
rd1 <- function(n) rnorm(n)
rd2 <- function(n) rnorMix(n,nM2)
rd3 <- function(n) rnorMix(n,nM3)
qd1 <- function(p) qnorm(p)
qd2 <- function(p) qnorMix(p,nM2)
qd3 <- function(p) qnorMix(p,nM3)
dd1 <- function(x) dnorm(x)
dd2 <- function(x) dnorMix(x,nM2)
dd3 <- function(x) dnorMix(x,nM3)

# Density plot to check D & T's Figure 1.
if(FALSE){
   x <- -500:500/100
   plot(x,dd1(x),type = "l", ylim = c(0,1))
   lines(x,dd2(x),col = 2)
   lines(x,dd3(x),col = 3)
   }

# Setup Monte-Carlo Experiment
#
require(doMC)
registerDoMC(3)

date()
system("hostname")
sessionInfo()
set.seed(1964)
opt <- list(set.seed = FALSE)


taus <- c(.1, .25,.5,.75,.9)
ns <- c(25,50,100,200)
dfs <- paste("d",1:3,sep = "")
R <- 200
A <- array(0,c(length(taus),2,length(ns),R))
B <- array(0,c(length(taus),4,length(ns),R))
AR <- array(0,c(length(ns),R))
FEJ <- foreach(j = 1:length(dfs),.options.multicore = opt) %dopar% {
      for(i in 1:length(ns)){
      makex <- parse(text = paste("r",dfs[j],"(",ns[i],")",sep = ""))
      for( k in 1:R){
	x <- eval(makex)
        f <- Qtile(x, taus)
	g <- QMCMC(x, taus, a0 = rep(0,5),S0 = diag(5))
	A[,1,i,k] <- f$Qhat 
	A[,2,i,k] <- g$Qhat 
	B[,,i,k] <- g$Cred
	AR[i,k] <- g$AcceptanceRatio
	}
     }
     list(A = A, B = B, AR = AR)
  }
