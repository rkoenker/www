# Neelon and Dunson Monotone Regression

rtnorm <- function(n,m,s,L,U) # Random truncated normals
    qnorm(runif(n,pnorm(L,m,s),pnorm(U,m,s)),m,s)

rtgmix <- function(n, m, s, delta){
   # Neelon and Dunson (Biometrics, 2004) truncated normal mixture
   # m,s are vectors corresponding to mean and sds of components
   # note that in the boundary cases there are only 2 components, not 3
   k <- length(m)
   vhat <- 1/sum(1/s^2)
   mhat <- vhat * sum(m/s^2)
   vtil <- 1/sum(1/s[-k]^2)
   mtil <- vtil * sum(m[-k]/s[-k]^2)
   A <- pnorm(delta,mtil,sqrt(vtil))/dnorm(0,mtil,sqrt(vtil))
   B <- (1- pnorm(delta,mhat,sqrt(vhat)))/dnorm(0,mhat,sqrt(vhat))
   prob <- A/(A+B)
   ifelse(runif(n) < prob, rtnorm(n,mtil,sqrt(vtil),-Inf,delta),
      rtnorm(n,mhat,sqrt(vhat),delta,Inf))
   }

PWLBasis <- function(x,g, intercept = TRUE){
   n <- length(x)
   p <- length(g)
   W <- matrix(0,n,p-1)
   for(i in 1:n){
      for(j in 2:p){
         W[i,j-1] <- max(min(x[i],g[j]) - g[j-1],0)
      }
   }
   if(intercept) W <- cbind(1,W)
   W
}
bmlm <- function(x,y,g = seq(0,1,length = 10), delta = .05, R = 1000, burn = 500){
   # Bayesian Monotone Bivariate Regression (Neelon and Dunson)
   n <- length(y)
   p <- length(g)
   # Design matrix for piecewise linear fit
   W <- PWLBasis(x,g, intercept = TRUE)
   # Prior for intercept parameter beta1~ N(mu1, sig1)
   mu1 <- 0
   sig1 <- sqrt(10)
   # Prior for initial slope parameter beta2~ N(mu2, sig2)
   mu2 <- 0
   sig2 <- sqrt(10)
   # Prior for precision parameter 1/sigma^2 ~ Gamma(nu1,nu2)
   nu1 <- nu2 <- 0.1
   # Prior for precision parameter lambda ~ Gamma(lam1, lam2)
   lam1 <- p/25
   lam2 <- 1
   # Prior for delta is unit point mass at 0.05 (provisionally)
   # Initialization:
   alpha <- beta <- rep(0,p)
   A <- matrix(0,p,R)
   # Main Loop:  
   for(i in 1:R){
      sigma <- sqrt(1/rgamma(1,nu1 + n/2, nu2 + crossprod(y - W %*% alpha)))
      lambda <- rgamma(1, lam1 + (p-1)/2, lam2 + sum(diff(beta[-1])^2)/2)
      # Intercept parameter called alpha by Neelon and Dunson
      psig1 <- sqrt(1/((1/sig1^2) + n/sigma^2))
      pmu1 <- (psig1^2) * (mu1/sig1^2 + sum(y - W[,-1] %*% alpha[-1])/sigma^2)
      alpha[1] <- beta[1] <- rnorm(1,pmu1,psig1)
      # Inital slope parameter
      bsig <- sqrt(c(sig2^2,  1/lambda, sigma^2/sum(W[,2]^2)))
      b2 <- lsfit(W[,2], y - W[,-2] %*% alpha[-2],intercept = FALSE)$coef
      bmu <- c(mu2, beta[3], b2)
      beta[2] <- rtgmix(1,bmu,bsig,delta)
      alpha[2] <- beta[2] * (beta[2] > delta)
      for(j in 3:p){
         bj <- lsfit(W[,j], y - W[,-j] %*% alpha[-j], intercept = FALSE)$coef
         if(j < p){
              bmu <- c(beta[j-1], beta[j+1], bj)
              bsig <- sqrt(c(1/lambda, 1/lambda, sigma^2/sum(W[,j]^2)))
              }
         else { #j = p, so there is no beta[j+1] to contribute
              bmu <- c(beta[j-1], bj)
              bsig <- sqrt(c(1/lambda, sigma^2/sum(W[,j]^2)))
              }
         beta[j] <- rtgmix(1,bmu,bsig,delta)
         alpha[j] <- beta[j] * (beta[j] > delta)
      }   
      A[,i] <- alpha
   }
   A[,(burn + 1):R]
}

x <- sort(runif(100))
y <- pbeta(x,3,3) + rnorm(100)/25
plot(x,y)
lines(x, pbeta(x,3,3), lwd = 2, )
W <- PWLBasis(x,g = seq(0,1,length = 10))
lines(x, lm(y ~ W - 1)$fitted, col = "red")
f <- bmlm(x,y)
lines(x,W %*% apply(f,1,mean),col = "blue")
