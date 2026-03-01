"rq.boot"<-
function(x, y, tau = 0.5, R = 1000, int = T, h = 0, z, method = "xy")
{
#bootstrap rq function --  returns matrix of R bootstrap estimates of rq vector
	n <- length(y)
	x <- as.matrix(x)
	if(int)
		X <- cbind(1, x)
	else X <- x
	p <- ncol(x) + int
	B <- matrix(0, R, p)
	if(method == "xy") {
#xy-pairs bootstrap
		s <- matrix(sample(n, n * R, replace = T), n, R)
		B <- xyboot(X, y, s, tau)
	}
	else if(method == "eqf") {
		if(missing(z))
			z <- rq(x, y, int = int)
		s <- sample(z$sol[2,  ], size = n * R, prob = c(diff(z$sol[1,  
			]), 0), replace = T)
		s <- matrix(s, n, R)
		bhat <- z$sol[ - (1:2), sum(z$sol[1,  ] < tau)]
		Y <- c(X %*% bhat) + s - mean(y)
		if(h > 0)
			Y <- Y + h * rnorm(n * R)
		B <- rqs(X, Y, tau)
	}
	else if(method == "heqf") {
		if(missing(z))
			z <- rq(x, y, int = int)
		m <- length(z$sol[1,  ]) - 1
		s <- matrix(sample(m, size = n * R, prob = diff(z$sol[1,  ]), 
			replace = T), n, R)
		bhat <- array(z$sol[ - (1:2), s], c(p, n, R))	
		Y <- heqf.Y(X, bhat)
		if(h > 0)
			Y <- Y + h * rnorm(n * R)
		B <- rqs(X, Y, tau)
	}
	else if(method== "resid") {
#residual bootstrap
		s <- matrix(sample(n, n * R, replace = T), n, R)
		if(missing(z))
			z <- rq(x, y, tau, int = int)
		Y <- c(X %*% z$coef) + matrix(z$resid[s], n, R)
		if(h > 0)
			Y <- Y + h * rnorm(n * R)
		B <- rqs(X, Y, tau)
	}
	else if(method == "pwy"){
		U <- t(X) %*% matrix((((runif(n * R) > tau) - tau)/n^0.5), n, R)
		B <- pwy(U,X,y,tau,R)
	}
	else
		stop("your chosen method is not allowed")
	B
}
"rqs"<-
function(x, y, theta = 0.5, tol = 0.0001)
{
#function to compute regression quantiles for multiple y's
#stripped down for monte-carlo purposes
	x <- as.matrix(x)
	p <- ncol(x)
	n <- nrow(x)
	k <- ncol(y)
	z <- .Fortran("rqs",
		as.integer(n),
		as.integer(p),
		as.integer(k),
		as.integer(n + 5),
		as.integer(p + 2),
		as.single(x),
		as.single(y),
		as.single(theta),
		as.single(tol),
		flag = as.integer(1),
		coef = single(p * k),
		resid = single(n),
		integer(n),
		single((n + 5) * (p + 2)),
		single(n),
		as.integer(1),
		as.integer(1),
		sol = single((p + 2)),
		dsol = single(n),
		lsol = as.integer(0),
		h = integer(p))
	return(t(matrix(z$coef, p, k)))
}
"xyboot"<-
function(x, y, s, theta = 0.5, tol = 0.0001)
{
#function to compute xypairs bootstrap for regression quantiles 
#stripped down for monte-carlo purposes
	x <- as.matrix(x)
	p <- ncol(x)
	n <- nrow(x)
	k <- ncol(s)
	z <- .Fortran("xys",
		as.integer(n),
		as.integer(p),
		as.integer(k),
		as.integer(n + 5),
		as.integer(p + 2),
		as.single(x),
		as.single(y),
		as.single(theta),
		as.single(tol),
		flag = as.integer(1),
		coef = single(p * k),
		resid = single(n),
		integer(n),
		single((n + 5) * (p + 2)),
		single(n),
		as.integer(1),
		as.integer(1),
		sol = single((p + 2)),
		dsol = single(n),
		lsol = as.integer(0),
		h = integer(p),
		xx = single(n * p),
		yy = single(n),
		as.integer(s))
	return(t(matrix(z$coef, p, k)))
}
"pwy"<-
function(U,X, y, tau = 0.5, N = 1000, tol=1e-4)
{
#resampling method of parzen,wei,ying for quantile regression
#NB. x should be full design matrix including intercept
	n <- length(y)
	p <- ncol(X)
	Y <- c(y,500000)
	x <- rbind(X,0)
	xu <- t(U)*n^.5/tau
	n <- n+1
	z<-.Fortran("pwy",
		as.integer(n),
		as.integer(p),
		as.integer(N),
		as.integer(n + 5),
		as.integer(p + 2),
		as.single(xu),
		as.single(x),
		as.single(Y),
		as.single(tau),
		as.single(tol),
		flag = as.integer(1),
		coef = single(p * N),
		resid = single(n),
		integer(n),
		single((n + 5) * (p + 2)),
		single(n),
		as.integer(1),
		as.integer(1),
		sol = single((p + 2)),
		dsol = single(n),
		lsol = as.integer(0),
		h=integer(p))
	return(t(matrix(z$coef, p, N)))
}
"heqf.Y"<-
function(x, b)
{
	n <- nrow(x)
	p <- ncol(x)
	R <- dim(b)[3]
	if(dim(b)[1] != p)
		stop("p dim don't match")
	if(dim(b)[2] != n)
		stop("n dim don't match")
	z <- .Fortran("heqfy",
		as.integer(n),
		as.integer(p),
		as.integer(R),
		as.single(x),
		as.single(b),
		y = single(n * R))
	matrix(z$y, n, R)
}
