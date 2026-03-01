rq.boot_function(x, y, tau = 0.5, int = F, h = 0, z = 0, method = "xy")
{
	n <- length(y)
	if(method == "xy") {
#xy-pairs bootstrap
		s <- sample(n, replace = T)
		rq(x[s,  ], y[s] + h * rnorm(n), tau, int = int)
	}
#eqf bootstrap
	if(method == "eqf") {
		s <- sample(z$sol[2,  ], size = n, prob = c(diff(z$sol[1,  ]),
			0), replace = T)
		bhat <- z$sol[ - (1:2), sum(z$sol[1,  ] < tau)]
		ys <- x %*% bhat + s - mean(y) + h * rnorm(n)
		rq(x, ys, tau, int = int)
	}
#heqf bootstrap
	if(method == "heqf") {
		m <- length(z$sol[1,  ]) - 1
		s <- sample(m, size = n, prob = diff(z$sol[1,  ]), replace = T)
			
		bhat <- z$sol[ - (1:2), s]
		ys <- diag(x %*% bhat) + h * rnorm(n)
		rq(x, ys, tau, int = int)
	}
	else {
#residual bootstrap
		s <- sample(n, replace = T)
		ys <- x %*% z$coef + z$resid[s] + h * rnorm(n)
		rq(x, ys, tau, int = int)
	}
}
