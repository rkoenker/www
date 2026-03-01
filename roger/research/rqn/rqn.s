"rqfn"<-
function(x, y, tau = 0.5, int = T, beta = 0.99995000000000001, eps = 
	1.0000000000000001e-05)
{
	if(!is.loaded(symbol.For("rqfn"))) dyn.load("src/rqfn/fn.o")	
	#rq function using frisch newton interior point algorithm
	n <- length(y)
	if(int)
		x <- cbind(1, x)
	x <- t(as.matrix(x))
	p <- nrow(x)
	if(n != ncol(x))
		stop("x and y don't match n")
	z <- .Fortran("rqfn",
		as.integer(n),
		as.integer(p),
		a = as.double(x),
		c = as.double( - y),
		rhs = double(p),
		d = double(n),
		beta = as.double(beta),
		eps = as.double(eps),
		tau = as.double(tau),
		wn = double(10 * n),
		wp = double((p + 3) * p),
		aa = double(p * p),
		it.count = integer(2))
	return(coef =  - z$wp[1:p], it.count = z$it.count)
}
"RQFN"<-
function(x, y, int = T, tau = 0.5, se = F, beta = 0.99995000000000001, eps = 
	1.0000000000000001e-05)
{
#new rq function for n large --RWK. 8.96. 
	if(!is.loaded(symbol.For("rqn"))) dyn.load("src/rqfn/fn.o")
	if(int)
		x <- cbind(1, x)
	x <- rbind(x, 0, 0)
	y <- c(y, 0, 0)
	p <- ncol(x)
	n2 <- length(y)
	mlim <- floor(10 * n2^(2/3))
	if(n2 != nrow(x))
		stop("x,y don't match n")
	z <- .Fortran("rqm",
		as.integer(n2),
		as.integer(p),
		as.double(t(x)),
		as.double( - y),
		double(p),
		as.double(rep(1, n2)),
		double(mlim * 10),
		wp = double(p * (p + 3)),
		as.double(beta),
		as.double(eps),
		as.double(tau),
		integer(n2),
		double(p * p),
		hist = integer(3 * 32))
	hist <- z$hist
	hist <- matrix(hist[hist > 0], 3)
	if(se) {
#quick and dirty se's
		n <- n2 - 2
		xxinv <- diag(p)
		xxinv <- backsolve(qr(x)$qr[1:p, 1:p], xxinv)
		xxinv <- xxinv %*% t(xxinv)
		b <-  - z$wp[1:p]
		r <- y - x %*% b
		pz <- sum(abs(r) < .Machine$single.eps)
		dn <- n * dn(tau, n, hs = F)
		ir <- (pz + 1):(dn + pz + 1)
		sparsity <- l1fit(ir/(n - p), sort(r[order(abs(r))][ir]))$coef[
			2]
		cov <- sparsity^2 * xxinv * tau * (1 - tau)
		return(coef = b, se = sqrt(diag(cov)))
	}
	else return(coef =  - z$wp[1:p], hist)
}
"monte"<-
function(run, ns = 1000 * c(40, 60, 80, 100, 120, 140, 160, 180), ps = c(4, 8), 
	R = 10, methods = expression(lm(y ~ x), rqfn(x, y), RQFN(x, y)), dfx = 
	expression(matrix(rnorm(p * n), n, p)), dfy = expression(rnorm(n)), mse
	 = list(2, 3))
{
	version <- 8	#function for timing experiments for rqn paper 
#Input:
#	ns-a vector of sample sizes
#	ps-a vector of parameter dimensions, intercept will be appended
#	R -number of replications of each n,p pair
#	methods-methods to be compared should be of the form:
#		expression(lm(y~x),l1fit(x,y),rqfn(x,y),RQFN(x,y))
#		this is a list which can be evaluated as eval(methods[[i]])
#	dfx-expression to generate design matrix
#	dfy-expression to generate response vectors
#	mse-list describing how to evaluate accuracy:
#		1.  benchmark method (number in methods list)
#		2.  new methods  under test (numbers in methods list)
#
#Output:
#	result-data structure with the components
#		times-array of timings
#		err  -root mse of bhat for new methods vis a vis benchmark
#		seed -initial .Random.seed
#	doc-attribute of result describing in detail how it was created	
#
	options(object.size = 150000000)	
	#checking for dynloading now occurs in the rq functions
#do the biggest problem first in case there are memory problems
	ns <- rev(sort(ns))
	ps <- rev(sort(ps))
	times <- array(0, c(length(methods), R, length(ps), length(ns)))	
	#err <- array(0, c(length(mse[[2]]), R, length(ps), length(ns)))
	seed <- .Random.seed
	for(i in 1:length(ns)) {
		n <- ns[i]
		print(paste("n=", n))
		for(j in 1:length(ps)) {
			p <- ps[j]
			print(paste("p=", p))
			b <- matrix(0, p + 1, length(methods))
			x <- eval(dfx)
			m <- x %*% rep(1, p)
			for(k in 1:R) {
				print(paste("k=", k))
				y <- m + eval(dfy)
				for(l in 1:length(methods)) {
				  times[l, k, j, i] <- unix.time(b[, l] <- eval(
				    methods[[l]])$coef)[1]
				}
#err[, k, j, i] <- sqrt(apply((b[, mse[[2]]] - b[
# , mse[[1]]])^2, 2, "mean"))
			}
		}
	}
	dimnames(times) <- list(paste(methods), NULL, paste("p=", ps, sep = ""),
		paste("n=", ns, sep = ""))
	result <- list(times = times, seed = seed, dfx = dfx, dfy = dfy)
	doc(result) <- how.created(paste("Test", run, "on", unix("hostname")), 
		text = F)
	return(result)
}
