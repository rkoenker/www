"rqx"<-
function(x, y, tau = 0.5, int = T, max.it = 50)
{
# Barrodale and Roberts quantile regression algorithm - lite
#			author:  R. Koenker, 
#                       date:    June, 1998
# The algorithm of Barrodale and Roberts (1974) for solving discrete l1 regression 
# problems has served as the standard simplex algorithm for quantile regression
# since the mid-70's.  For large problems recent developments of interior point
# methods have proven to be superior, but for small to moderate problems the
# approach of Barrodale and Roberts is still very competitive.  Unfortunately,
# the fortran version of original algorithm is somewhat obscure.  In this pure
# R version I have tried to reduce the algorithm to its essentials -- at each
# iteration we find the Cauchy direction i.e. the direction of steepest descent
# and then we take this direction until it no longer reduces the value of the
# of the objective function.  This "Cauchy Step" involves solving a one dimensional
# "through the origin" rq problem.  This can be formulated as a weighted quantile
# problem that is solved by the function wquantile().  A more efficient version
# of this algorithm should eventually be coded in fortran with an O(n) algorithm
# for solving these subproblems using the approach of Floyd and Rivest (1975), if
# possible.  Meanwhile this is simply a heuristic device.
#
	if(int) x <- cbind(1, x)
	p <- ncol(x)
	n <- nrow(x)	#Phase I -- find a random (!) initial basis
	h <- sample(1:n, size = p)
	it <- 0
	repeat {
		it <- it + 1
		Xhinv <- solve(x[h,  ])
		bh <- Xhinv %*% y[h]
		rh <- y - x %*% bh	
	#find direction of steepest descent along one of the edges
		g <-  - t(Xhinv) %*% t(x[ - h,  ]) %*% c(tau - (rh[ - h] < 0))
		g <- c(g + (1 - tau),  - g + tau)
		ming <- min(g)
		if(ming >= 0 || it > max.it)
			break
		h.out <- seq(along = g)[g == ming]
		sigma <- ifelse(h.out <= p, 1, -1)
		if(sigma < 0)
			h.out <- h.out - p
		d <- sigma * Xhinv[, h.out]	
	#find step length by one-dimensional minimization
		xh <- x %*% d
		step <- wquantile(xh, rh, tau)
		h.in <- step$k
		h <- c(h[ - h.out], h.in)
	}
	if(it > max.it)
		warning("non-optimal solution: max.it exceeded")
	return(bh)
}
"rho.tau"<- function(r, tau = 0.5) tau * pmax(r, 0) + (1 - tau) * pmax( - r, 0)
"wquantile"<- function(x, y, t = 0.5) {
#weighted quantile by brute force (full sorting) for scalar quantile regression
# through the origin model:  Q_{y_i} (t)  = x_i b
#
	ord <- order(y/x)
	b <- (y/x)[ord]
	wabs <- abs(x[ord])
	k <- sum(cumsum(wabs) < ((t - 0.5) * sum(x) + 0.5 * sum(wabs)))
	list(b = b[k + 1], k = ord[k + 1])
}
