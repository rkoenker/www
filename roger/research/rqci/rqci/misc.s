"as.single.bug"<-
function()
{
	.Random.seed <- c(9, 60, 58, 4, 5, 1, 62, 31, 56, 63, 11, 2)
	x <- cbind(rnorm(50), rbinom(50, 1, 0.5), runif(50))
	y <- rnorm(50)
	.Random.seed <- c(41, 59, 49, 42, 57, 2, 40, 55, 31, 32, 6, 1)
	n <- 50
	R <- 1000
	tau <- 0.5
	X <- cbind(1, x)
	U <- t(X) %*% matrix((((runif(n * R) > tau) - tau)/n^0.5), n, R)
	xu <- (t(U) * n^0.5)/tau
	return(xu[3591])
}
"boot.ci"<-
function(B, alpha = 0.10000000000000001)
{
#compute percentile confidence intervals from bootstrap sample matrix B
	apply(B, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
}
"ci.test"<-
function(n)
{
	rq.ci(c(rep(-1, 2 * n + 1), rep(1, 2 * n + 1)), rep( - n:n, 2))
}
"d0"<-
function(p, n, hs = F)
{
	PI <- 3.1415899999999999
	x0 <- qnorm(p)
	f0 <- (1/sqrt(2 * PI)) * exp(.Uminus((x0^2/2)))
	if(hs == T)
		n^(-1/3) * qnorm(1 - 0.025000000000000001)^(2/3) * ((1.5 * f0^2
			)/(2 * x0^2 + 1))^(1/3)
	else n^-0.20000000000000001 * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^
			0.20000000000000001
}
"dn"<-
function(p, n, alpha = 0.10000000000000001, hs = F)
{
	PI <- 3.1415899999999999
	x0 <- qnorm(p)
	f0 <- (1/sqrt(2 * PI)) * exp(.Uminus((x0^2/2)))
	if(hs == T)
		dn <- n^(-1/3) * qnorm(1 - alpha/2)^(2/3) * ((1.5 * f0^2)/(2 * 
			x0^2 + 1))^(1/3)
	else dn <- n^-0.20000000000000001 * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^
			0.20000000000000001
	if(min(p - dn) < 0 | max(p + dn) > 1) {
		warning("n too small for specified p")
		return((p <= 0.5) * pmin(p, dn) + (p > 0.5) * pmin(1 - p, dn))
	}
	else return(dn)
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
"monte"<-
function(n = 50, p = 3, dfx = 3, dfy = 3, N = 1000, memory = F, methods = c(
	"score", "sparsity", "boot.pwy", "boot.heqf", "boot.xy"))
{
#Monte-carlo function for Prague rqci paper
#n=sample size
#p=dimension of model regression parameter minus one for intercept
#dfx= df's of the t distrution used to generate the x's (iid-ly) -- fixed
#dfx= df's of the t distrution used to generate the y's (iid-ly)
#N= number of replications
#Returns:
#	seed = random seed needed to replicate experiment
#	A = Array of results (2 by p+1 by 3 by N) with contents as follows:
#		upper and lower limits
#		coeficients
#		version:  
#			1.  ranktestinversion, 
#			2.  sparsity-hall-sheather, 
#			3.  bootstrap pwy
#			4.  bootstrap heqf
#			5.  bootstrap xy
#		replications		
#Use sum.A to summarize the results in component A.
	seed <- .Random.seed
	if(memory == T)
		mem <- rep(0, N)
	x <- matrix(rt(n * p, dfx), n, p)
	qn <- 1/diag(solve(crossprod(cbind(1, x))))
	y <- matrix(rt(n * N, dfy), n, N)
	A <- array(0, dim = c(2, p + 1, length(methods), N))
	for(i in 1:N) {
		for(j in 1:length(methods)) {
			method <- methods[j]
			A[,  , j, i] <- rq.ci(x, y[, i], qn = qn, method = 
				method)
			if(memory == T)
				mem[i] <- memory.size()
		}
	}
	if(memory == T)
		return(seed, A, mem)
	else return(seed, A)
}
"monte.het"<-
function(n = 50, p = 3, dfy = 3, N = 1000, memory = F, methods = c("score", 
	"sparsity", "boot.pwy", "boot.heqf", "boot.xy"))
{
#Monte-carlo function for Prague rqci paper heteroscedastic case
#n=sample size
#p=dimension of model regression parameter minus one for intercept
#dfx= df's of the t distrution used to generate the y's (iid-ly)
#N= number of replications
#Returns:
#	seed = random seed needed to replicate experiment
#	A = Array of results (2 by p+1 by 3 by N) with contents as follows:
#		upper and lower limits
#		coeficients
#		version:  
#			1.  ranktestinversion, 
#			2.  sparsity-hall-sheather, 
#			3.  bootstrap pwy
#			4.  bootstrap heqf
#			5.  bootstrap xy
#		replications		
#Use sum.A to summarize the results in component A.
	seed <- .Random.seed
	if(memory == T)
		mem <- rep(0, N)
	x <- matrix(exp(rnorm(n * p)), n, p)
	qn <- 1/diag(solve(crossprod(cbind(1, x))))
	scale <- (1 + apply(x, 1, "sum"))/5
	y <- matrix(rt(n * N, dfy), n, N) * scale
	A <- array(0, dim = c(2, p + 1, length(methods), N))
	for(i in 1:N) {
		for(j in 1:length(methods)) {
			method <- methods[j]
			A[,  , j, i] <- rq.ci(x, y[, i], qn = qn, method = 
				method)
			if(memory == T)
				mem[i] <- memory.size()
		}
	}
	if(memory == T)
		return(seed, A, mem)
	else return(seed, A)
}
"os"<-
function(mode = "function", where = 1, long = F)
{
	z <- objects.summary(where = where, mode = "function")
	if(long)
		return(z)
	else return(row.names(z))
}
"pwy"<-
function(U, X, y, tau = 0.5, N = 1000, tol = 0.0001)
{
#resampling method of parzen,wei,ying for quantile regression
#NB. x should be full design matrix including intercept
	n <- length(y)
	p <- ncol(X)
	Y <- c(y, 500000)
	x <- rbind(X, 0)
	xu <- (t(U) * n^0.5)/tau
	n <- n + 1
	z <- .Fortran("pwy",
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
		h = integer(p))
	return(t(matrix(z$coef, p, N)))
}
"pwyr"<-
function(X, y, tau = 0.5, N = 1000, tol = 0.0001)
{
#resampling method of parzen,wei,ying for quantile regression
#NB. x should be full design matrix including intercept
	n <- length(y)
	p <- ncol(X)
	U <- t(X) %*% matrix((((runif(n * N) > tau) - tau)/n^0.5), n, N)
	Y <- c(y, 500000)
	x <- rbind(X, 0)
	xu <- rbind(X, (t(U) * n^0.5)/tau)
	n <- n + 1
	z <- .Fortran("pwy",
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
		lsol = as.integer(0))
	return(t(matrix(z$coef, p, N)))
}
"rq.boot"<-
function(x, y, tau = 0.5, alpha = 0.10000000000000001, R = 500, int = T, h = 0, 
	z, method = "boot.xy")
{
#bootstrap rq function --  returns matrix of R bootstrap estimates of rq vector
	n <- length(y)
	x <- as.matrix(x)
	if(int)
		X <- cbind(1, x)
	else X <- x
	p <- ncol(x) + int
	B <- matrix(0, R, p)
	if(method == "boot.xy") {
#xy-pairs bootstrap
		s <- matrix(sample(n, n * R, replace = T), n, R)
		B <- xyboot(X, y, s, tau)
	}
	else if(method == "boot.eqf") {
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
	else if(method == "boot.heqf") {
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
	else if(method == "boot.resid") {
#residual bootstrap
		s <- matrix(sample(n, n * R, replace = T), n, R)
		if(missing(z))
			z <- rq(x, y, tau, int = int)
		Y <- c(X %*% z$coef) + matrix(z$resid[s], n, R)
		if(h > 0)
			Y <- Y + h * rnorm(n * R)
		B <- rqs(X, Y, tau)
	}
	else if(method == "boot.pwy") {
		U <- t(X) %*% matrix((((runif(n * R) > tau) - tau)/n^0.5), n, R
			)
		B <- pwy(U, X, y, tau, R)
	}
	else stop("your chosen method is not allowed")
	apply(B, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
}
"rq.boot.xy"<-
function(x, y, tau = 0.5, R = 5, int = T, h = 0, z = 0, method = "xy")
{
#bootstrap rq function --  returns matrix of R bootstrap estimates of rq vector
	n <- length(y)
	x <- as.matrix(x)
	p <- ncol(x) + int
	B <- matrix(0, R, p)
	if(method == "xy") {
		browser()	#xy-pairs bootstrap
		s <- matrix(sample(n, n * R, replace = T), n, R)
		Y <- matrix(y[s], n, R)
		if(h > 0)
			Y <- Y + h * matrix(rnorm(n * R), n, R)
		for(i in 1:R)
			B[i,  ] <- rq(x[s[, i],  ], Y[, i], tau, int = int)$
				coef
	}
	B
}
"rq.ci"<-
function(x, y, theta = 0.5, qn, alpha = 0.10000000000000001, int = T, tol = 
	0.0001, method = "score", hs = T, interpolate = T, tcrit = T)
{
#function to compute regression quantiles
	if(int) {
		x <- cbind(1, x)
	}
	big <- .Machine$single.xmax
	x <- as.matrix(x)
	p <- ncol(x)
	n <- nrow(x)
	if(n != length(y))
		stop("x and y don't match n")
	if(tcrit)
		cutoff <- qt(1 - alpha/2, n - p)
	else cutoff <- qnorm(1 - alpha/2)
	if(missing(qn))
		qn <- 1/diag(solve(crossprod(x)))
	z <- rqci(x, y, theta, qn, cutoff)
	if(method == "score") {
		Tci <- z$ci
		Tn <- z$Tn
		Tci[3,  ] <- Tci[3,  ] + (cutoff - abs(Tn[3,  ]))/abs(Tn[4,  ] - 
			Tn[3,  ]) * abs(Tci[4,  ] - Tci[3,  ])
		Tci[2,  ] <- Tci[2,  ] - (cutoff - abs(Tn[2,  ]))/abs(Tn[1,  ] - 
			Tn[2,  ]) * abs(Tci[1,  ] - Tci[2,  ])
		Tci[2,  ][is.na(Tci[2,  ])] <-  - big
		Tci[3,  ][is.na(Tci[3,  ])] <- big
		return(ci = Tci[2:3,  ])
	}
	else if(method == "sparsity") {
		all.rqs <- rq(x, y, int = F)
		isol <- sum(all.rqs$sol[1,  ] <= theta)
		Tcoef <- all.rqs$sol[3:(p + 2), isol]
		se.hs <- sqrt(rq.omega(all.rqs, theta, alpha, hs = T)/qn)
		ci.hs <- rbind(Tcoef - cutoff * se.hs, Tcoef + cutoff * se.hs)
		return(ci.hs)
	}
	else if(charmatch("boot", method, nomatch = F) != F) {
		if(!exists("all.rqs", frame = 0) & length(grep("eqf", method)) > 
			0)
			all.rqs <- rq(x, y, int = F)
		return(rq.boot(x, y, theta, alpha, z = all.rqs, int = F, method
			 = method))
	}
}
"rq.omega"<-
function(z, tau = 0.5, alpha = 0.10000000000000001, hs = F, tol = 0.0001)
{
#computes omega^2=tau(1-tau)/f^2(xi_a) for se of rq's
#uses linearized eqf from the full rq-process see function qrq()
#bandwidth is either Bofinger (hs=F) or Hall-Sheather (hs=T)
#z is the rq structure 
	dn <- dn(tau, nrow(z$dsol), alpha = alpha, hs = hs)
	q <- qrq(z, c(tau - dn, tau + dn))
	shat <- diff(q)/(2 * dn)
	tau * (1 - tau) * shat^2
}
"rqci"<-
function(x, y, theta, qn, cutoff, tol = 0.0001)
{
#function to compute regression quantiles
	big <- .Machine$single.xmax
	p <- ncol(x)
	n <- nrow(x)
	z <- .Fortran("rqci",
		as.integer(n),
		as.integer(p),
		as.integer(n + 5),
		as.integer(p + 2),
		as.integer(p + 4),
		as.single(x),
		as.single(y),
		as.single(theta),
		as.single(tol),
		flag = as.integer(1),
		coef = single(p),
		resid = single(n),
		integer(n),
		single((n + 5) * (p + 4)),
		single(n),
		as.integer(2),
		as.integer(2),
		sol = single(2 * p),
		dsol = single(2 * n),
		lsol = as.integer(0),
		h = integer(p),
		as.single(qn),
		as.single(cutoff),
		ci = single(4 * p),
		tnmat = single(4 * p),
		as.single(big),
		as.logical(T))
	return(coef = z$coef, ci = matrix(z$ci, 4), Tn = matrix(z$tnmat, 4))
}
"rqci.new"<-
function(x, y, theta = 0.5, alpha = 0.10000000000000001, int = T, tol = 0.0001, 
	method = "inv.score", interpolate = T, tcrit = T)
{
#function to compute regression quantiles
	if(int) {
		x <- cbind(1, x)
	}
	big <- .Machine$single.xmax
	x <- as.matrix(x)
	p <- ncol(x)
	n <- nrow(x)
	if(n != length(y))
		stop("x and y don't match n")
	if(tcrit)
		cutoff <- qt(1 - alpha/2, n - p)
	else cutoff <- qnorm(1 - alpha/2)
	qn <- 1/diag(solve(crossprod(x)))
	z <- .Fortran("rqci",
		as.integer(n),
		as.integer(p),
		as.integer(n + 5),
		as.integer(p + 2),
		as.integer(p + 4),
		as.single(x),
		as.single(y),
		as.single(theta),
		as.single(tol),
		flag = as.integer(1),
		coef = single(p),
		resid = single(n),
		integer(n),
		single((n + 5) * (p + 4)),
		single(n),
		as.integer(nsol),
		as.integer(ndsol),
		sol = single((p + 2) * nsol),
		dsol = single(n * ndsol),
		lsol = as.integer(0),
		h = integer(p * nsol),
		qn = as.single(qn),
		cutoff = as.single(cutoff),
		ci = single(4 * p),
		tnmat = single(4 * p),
		as.single(big),
		as.logical(T))
	if(z$flag != 0)
		warning(switch(z$flag,
			"Solution may be nonunique",
			"Premature end - possible conditioning problem in x."))
	if(ci & interpolate & ci.method == "inv.score") {
		Tn <- matrix(z$tnmat, nrow = 4)
		Tci <- matrix(z$ci, nrow = 4)
		Tci[3,  ] <- Tci[3,  ] + (cutoff - abs(Tn[3,  ]))/abs(Tn[4,  ] - 
			Tn[3,  ]) * abs(Tci[4,  ] - Tci[3,  ])
		Tci[2,  ] <- Tci[2,  ] - (cutoff - abs(Tn[2,  ]))/abs(Tn[1,  ] - 
			Tn[2,  ]) * abs(Tci[1,  ] - Tci[2,  ])
		Tci[2,  ][is.na(Tci[2,  ])] <-  - big
		Tci[3,  ][is.na(Tci[3,  ])] <- big
		dimnames(Tci) <- list(c("Lower Bound", "Lower Bound", 
			"Upper Bound", "Upper Bound"), NULL)
		return(coef = z$coef, resid = z$resid, h = z$h[1:p], ci = Tci[2:
			3,  ])
	}
	if(ci & ci.method == "inv.score") {
		Tci <- matrix(z$ci, nrow = 4)
		dimnames(Tci) <- list(c("Lower Bound", "Lower Bound", 
			"Upper Bound", "Upper Bound"), NULL)
		return(coef = z$coef, ci = Tci, Tn = matrix(z$tnmat, nrow = 4))
	}
	else if(ci.method == "band") {
		z <- rq(x, y, int = F)
		se <- sqrt(rq.omega(z, theta)/qn)
		Tci <- matrix(0, 2, p)
		dimnames(Tci) <- list(c("Lower Bound", "Upper Bound"), NULL)
		isol <- sum(z$sol[1,  ] <= theta)
		Tcoef <- z$sol[3:(p + 2), isol]
		Tci[1,  ] <- Tcoef - cutoff * se
		Tci[2,  ] <- Tcoef + cutoff * se
		return(coef = Tcoef, ci = Tci)
	}
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
"sum.A"<-
function(a, trim = 0.050000000000000003)
{
#summarize output from monte()
#cov is the coverage frequency
#cil is mean confidence interval length
#NB only the slope parameters are used.
	a <- a$A
	p <- dim(a)[2] - 1
	m <- dim(a)[3]
	n <- dim(a)[4]
	cov <- rep(0, m)
	cil <- rep(0, m)
	for(i in 1:m) {
		cov[i] <- sum(a[1, -1, i,  ] * a[2, -1, i,  ] > 0)/(p * n)
		cil[i] <- mean(a[2, -1, i,  ] - a[1, -1, i,  ], trim = trim)
	}
	return(cov, cil)
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
