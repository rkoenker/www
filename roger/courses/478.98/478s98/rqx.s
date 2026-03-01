"rqx"<-
function(x, y, tau = 0.5, int = T, max.it = 50)
{
#Barrodale and Roberts quantile regression algorithm - lite
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
	#find direction of steepest descent along one the edges
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
		step <- wquantile(rh, xh, tau)
		h.in <- step$k
		h <- c(h[ - h.out], h.in)
	}
	if(it > max.it)
		warning("non-optimal solution: max.it exceeded")
	return(bh)
}
"wquantile"<-
function(x, y, t = 0.5)
{
#weighted quantile by brute force for model Q_{y_i} (t)  = x_i b
#
	ord <- order(y/x)
	b <- (y/x)[ord]
	wabs <- abs(x[ord])
	k <- sum(cumsum(wabs) < ((t - 0.5) * sum(w) + 0.5 * sum(wabs)))
	return(b = b[k + 1], k = ord[k + 1])
}
