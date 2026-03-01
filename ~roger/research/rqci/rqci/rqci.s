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
		return(coef = z$coef, ci = Tci[2:3,  ])
	}
	else if(method == "sparsity") {
		z <- rq(x, y, int = F)
		se <- sqrt(rq.omega(z, theta, hs = hs)/qn)
		Tci <- matrix(0, 2, p)
		isol <- sum(z$sol[1,  ] <= theta)
		Tcoef <- z$sol[3:(p + 2), isol]
		return(coef = Tcoef, ci = rbind(Tcoef - cutoff * se, Tcoef + 
			cutoff * se))
	}
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
		as.integer(1),
		as.integer(1),
		sol = single(1),
		dsol = single(1),
		lsol = as.integer(0),
		h = integer(1),
		as.single(qn),
		as.single(cutoff),
		ci = single(4 * p),
		tnmat = single(4 * p),
		as.single(big),
		as.logical(T))
	return(coef = z$coef, ci = matrix(z$ci, 4), Tn = matrix(z$tnmat, 4))
}
"rq.omega"<-
function(z, tau = 0.5, hs = F, tol = 0.0001)
{
#computes omega^2=tau(1-tau)/f^2(xi_a) for se of rq's
#uses linearized eqf from the full rq-process see function qrq()
#bandwidth is either Bofinger (hs=F) or Hall-Sheather (hs=T)
#z is the rq structure 
	dn <- dn(tau, nrow(z$dsol), hs = hs)
	q <- qrq(z, c(tau - dn, tau + dn))
	shat <- diff(q)/(2 * dn)
	tau * (1 - tau) * shat^2
}
"dn"<-
function(p, n, hs = F)
{
	PI <- 3.1415899999999999
	x0 <- qnorm(p)
	f0 <- (1/sqrt(2 * PI)) * exp(.Uminus((x0^2/2)))
	if(hs == T)
		dn <- n^(-1/3) * qnorm(1 - 0.025000000000000001)^(2/3) * ((1.5 * 
			f0^2)/(2 * x0^2 + 1))^(1/3)
	else dn <- n^-0.20000000000000001 * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^
			0.20000000000000001
	if(min(p - dn) < 0 | max(p + dn) > 1) {
		warning("n too small for specified p")
		return((p <= 0.5) * pmin(p, dn) + (p > 0.5) * pmin(1 - p, dn))
	}
	else return(dn)
}
"qrq"<-
function(s, a)
{
#computes linearized quantiles from rq data structure
#v is the rq structure e.g. rq(x,y)
#a is a vector of quantiles required
	if(min(a) < 0 | max(a) > 1) stop("alphas out of range [0,1]")
	r <- s$sol[1,  ]
	q <- s$sol[2,  ]
	q <- c(q[1], q)
	J <- length(r)
	r <- c(0, (r[1:J - 1] + r[2:J])/2, 1)
	u <- rep(0, length(a))
	for(k in 1:length(a)) {
		i <- sum(r < a[k])
		w <- (a[k] - r[i])/(r[i + 1] - r[i])
		u[k] <- w * q[i + 1] + (1 - w) * q[i]
	}
	u
}
