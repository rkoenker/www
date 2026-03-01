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
