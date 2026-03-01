# Example 1 from Candes:  Highly Robust Error Correction by Convex Programming
# Revised Version:  June 23, 2007 RWK.

DSelector <- function(X,y,K, sparse = FALSE){
# This is a simple implementation of the Dantzig selector procedure of Candes-Tao. 
# The basic idea of the Dantzig selector is to solve problems that look like this:
#
#	min || b ||_1 such that || A'r ||_oo < K
#
# where r = y - Ab.  But in this version we've implement a slightly more general
#version for:
#
#	min || b ||_1 such that -K <=  A'r <=  K
#
# for vector valued K.  We can formulate this in the usual QR setup as
#
#	rq.fit.fnc(X,y,R=R, r=r)
#
# which solves problems of the form:
#
# 	min || y - Xb ||_1 s.t. Rb >= r
#
# If X is sparse it would be good to use the sparse version, but this needs work.
#
require(quantreg)
n <- nrow(X)
p <- ncol(X)

A <- crossprod(X)
R <- rbind(A, -A)
a <- crossprod(X,y)
r <- c(a-K, -a-K)
zp <- rep(0,p)
if(sparse){
	Ip <- as(p,"matrix.diag.csr")
	f <- rq.fit.sfnc(Ip, zp, R=R, r=r)
	}
else{
	Ip <- diag(p)
	f <- rq.fit.fnc(Ip, zp, R=R, r=r)
	}
return(f)
}


LPDecoderRing <- function(A,y,sigma,Kval = 3){
# R Implementation of the Candes-Gauss-Dantzig Decoder Ring

   NullSpace <- function(A){
        S <- svd(A,nrow(A))
        r <- length(S$d)
        S$u[,(r+1):nrow(A),drop=FALSE]
       }
   Q <- t(NullSpace(A))
   Qy <- Q %*% y
   K <- sqrt(diag(crossprod(Q))) * sigma * Kval
   f <- DSelector(Q, Qy, K) 

# Reprojection Step

   s <- abs(f$coef) > sigma
   e <- rep(0,nrow(A))
   e[s] <- lm(Qy ~ Q[,s] -1)$coef
   xhat <- lm((y-e) ~ A - 1)$coef
   return(xhat)
  }

# Test Problem as in Section 4 of Candes (2007)
n <- 256
m <- 2 * n
A <- matrix(rnorm(m*n),m,n)
rho <- 0.10
k <- round(rho * m)
x <- rnorm(n)
y <- A %*% x
yt <- y
sigma <- median(abs(y))/16
s <- sample(1:m,k)
yt[s] <- -yt[s] # Original Candes example
u <- rnorm(m)*sigma
yt <- yt + u

xhat <- LPDecoderRing(A, yt, sigma)

pdf("Eg3.pdf",width=8,height=3)
par(mfrow=c(1,3))
g <- lm(yt ~ A -1)$coef
rmseG <- round(sqrt(var(x-g)),3)
plot(x,g,ylab=expression(hat(x)),main="Gauss")
text(1,-2,paste("RMSE = ",rmseG))
plot(x,xhat,ylab=expression(hat(x)),main="Dantzig")
rmseD <- round(sqrt(var(x-xhat)),3)     
text(1,-2,paste("RMSE = ",rmseD))
h <- lm((y + u) ~ A - 1)$coef
plot(x,h,ylab=expression(hat(x)),main="Oracle")
rmseH <- round(sqrt(var(x-h)),3)        
text(1,-2,paste("RMSE = ",rmseH))
dev.off()



