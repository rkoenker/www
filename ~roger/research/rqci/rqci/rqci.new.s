function(x, y, theta = 0.5, alpha = 0.1, int = T, tol = 0.0001, method = 
	"inv.score", interpolate = T, tcrit = T)
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
