DSelector <- function(X,y,sigma, lambda = 3.5, sparse = FALSE){
#This is a simple implementation of the Dantzig selector procedure of Candes-Tao, see
#the preprint pdf in the pdfs directory.  The Dantzig selector solves the problems 
#
#	min || b ||_1 such that || X'r ||_oo < K
#
#where K = lambda * sigma with sigma the scale of the noise and
#r = y - Xb.  So we can formulate this in the usual QR setup (with default tau = .5)
#
#	rq.fit.fnc(X,y,R=R, r=r)
#	rq.fit.sfnc(X,y,R=R, r=r)
#
#where y = rep(0,p), X = diag(p), A = X'X, R= rbind(A,-A), r = c(a - K, -a -K)
#and a = X'y.  In the latter version it is assumed that X is stored in 
#SparseM (matrix.csr) form.
#
#
require(quantreg)
n <- nrow(X)
p <- ncol(X)

K <- lambda * sigma
A <- t(X) %*% X
R <- rbind(A, -A)
a <- c(as.matrix(t(X) %*% y))
r <- c(a-K, -a-K)
zp <- rep(0,p)
if(sparse){
	Ip <- as(p,"matrix.diag.csr")
	f <- rq.fit.sfnc(Ip, zp, R=R, r=r)
	}
else{
	Ip <- diag(p)
	f <- rqc(Ip, zp, R=R, r=r)
	}
return(f)
}
"rqc" <- function(X,y,R,r,tau=.5){
    n <- nrow(X)
    p <- nrow(X)
    m <- nrow(R)
    u <- rep(1,n)
    a1 <- (1-tau)*u
    a2 <- rep(1,m)
    b <- t(X) %*% a1
    f <- lpfnc(t(X), -y, t(R), -r, b, u, a1, a2)
    return(coef = -f$coef, it = f$it)
}

