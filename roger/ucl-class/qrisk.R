"qrisk" <-
function(x, alpha=c(.1,.3), w = c(.7,.3), mu = .07, R = NULL, r = NULL, lambda = 10000){
#
# find optimal Choquet-risk portfolios given:
#
#	x 	(n by p) matrix of asset returns
#	alphas	alphas defining a Choquet capacity risk function
#	w	w defining weights for Choquet capacity risk function
#	R	Matrix defining constraints on the parameters 
#	r	rhs defining constraints on the parameters 
#	mu	required mean rate of return	
#	lambda	Lagrange multiplier for RoR constraint 
#
	n <- nrow(x)
	p <- ncol(x)
	m <- length(alpha)
	if(length(w)!=m)return("length of w doesn't match length of alpha")
	xbar <- apply(x,2,mean)
	y <- x[,1]
	r <- c(r,lambda*(xbar[1]-mu), -lambda*(xbar[1]-mu))
	X <- x[,1]-x[,-1]
	R <- rbind(R,lambda*(xbar[1]-xbar[-1]), -lambda*(xbar[1]-xbar[-1]))
	R <- cbind(matrix(0,nrow(R),m),R)
	fit <- rq.fit.hogg(X,y,taus=alpha,w=w,R=R,r=r)
	pihat <- c(1-sum(fit[-(1:m)]),fit[-(1:m)])
	yhat <- x%*%pihat
	etahat <- quantile(yhat,alpha)
	muhat <- mean(yhat)
	qrisk <- 0
	for(i in 1:length(alpha))
		qrisk <- qrisk + w[i]*sum(yhat[yhat<etahat[i]])/(n*alpha[i]) 
	return(pihat,muhat,qrisk)
}
"srisk" <-
function(x,mu=.07,lambda=100000000,alpha=.1,eps=.0001){
#
# find optimal sigma-risk (Markowitz) portfolios given:
#
#	x 	(n by p) matrix of asset returns
#	mu	required mean return
#	lambda	Lagrange multiplier for required mean return constraint
#
	n <- nrow(x)
	p <- ncol(x)
	X <- rbind(x,apply(x,2,mean))
	y <- X[,1]
	y[n+1] <- lambda*(y[n+1]-mu)
	X <- X[,1]-X[,-1]
	X <- cbind(c(rep(1,n),0),X)
	X[n+1,] <- lambda*X[n+1,]
	fit <- lm(y~X-1)
	pihat <- c(1-sum(fit$coef[-1]),fit$coef[-1])
	if(abs(fit$residual[n+1]) > eps) return("lambda too small?")
	yhat <- x%*%pihat
	muhat <- mean(x%*%pihat)
	sigma <- sqrt(var(x%*%pihat))
	return(pihat,muhat,sigma)
}
"rq.fit.hogg" <-
function (x, y, taus = c(.1,.3,.5), weights = c(.7,.2,.1),  
	R= NULL, r = NULL, beta = 0.99995, eps = 1e-06) 
{
    n1 <- length(y)
    n2 <- nrow(R)
    m <- length(taus)
    p <- ncol(x)+m
    if (n != nrow(x)) 
        stop("x and y don't match n")
    if (m != length(weights)) 
        stop("taus and weights differ in length")
    if (any(taus < eps) || any(taus > 1 - eps)) 
        stop("taus outside (0,1)")
    if (any(weights < eps) || any(weights > 1 - eps)) 
        stop("weights outside (0,1)")
    X <- cbind(kronecker(diag(m),rep(1,n)),kronecker(rep(1,m),x))
    y <- rep(y,m)
    rhs <- c(weights*(1 - taus)*n, sum(weights*(1-taus)) * apply(x, 2, sum))
    if(n2!=length(r))
	stop("R and r of incompatible dimension")
    if(ncol(R)!=p)
	stop("R and X of incompatible dimension")
    d <- rep(1, m*n)
    u <- rep(1, m*n)
    if(length(r)){
       wn1 <- rep(0, 10 * m*n)
       wn1[1:(m*n)] <- .5
       wn2 <- rep(0,6*n2)
       wn2[1:n2] <- 1 
       z <- .Fortran("rqfnc", as.integer(m*n), as.integer(n2), as.integer(p), 
           a1 = as.double(t(as.matrix(X))), c1 = as.double(-y), 
           a2 = as.double(t(as.matrix(R))), c2 = as.double(-r), 
           rhs = as.double(rhs), d1 = double(m*n), d2 = double(n2), 
           as.double(u), beta = as.double(beta), eps = as.double(eps), 
           wn1 = as.double(wn1), wn2 = as.double(wn2), wp = double((p + 3) * p), 
	   it.count = integer(2), info = integer(1), PACKAGE = "quantreg")
	}
    else{
	wn <- rep(0, 10 * m*n)
    	wn[1:(m*n)] <- .5
    	z <- .Fortran("rqfnb", as.integer(m*n), as.integer(p), a = as.double(t(as.matrix(X))), 
		c = as.double(-y), rhs = as.double(rhs), d = as.double(d), as.double(u), 
		beta = as.double(beta), eps = as.double(eps), wn = as.double(wn), 
		wp = double((p + 3) * p), it.count = integer(2), info = integer(1))
	}
    if (z$info != 0) 
        stop(paste("Error info = ", z$info, "in stepy: singular design"))
    coefficients <- -z$wp[1:p]
    return(coefficients)
}
