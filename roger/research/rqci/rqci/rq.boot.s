"rq.boot"<-
function(x, y, tau = 0.5, R = 1000, int = T, h = 0, z = 0, method = "xy")
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
		B <- xyboot(X, y, s)
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
	#Y <- matrix(X[i, j] %tensor% bhat[j, k, l])
#Y <- matrix(Y[rep(diag(n), R) == 1], n, R)
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
		B <- pwy(X,y,tau,R)
	}
	else
		stop("your chosen method is not allowed")
	B
}
