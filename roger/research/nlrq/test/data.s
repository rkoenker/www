"bard"<-
function(t)
{
#bard model for nonlinear l1 estimation
	x <- x.bard
	y <- y.bard
	f <- y - t[1] - x[, 1]/(t[2] * x[, 2] + t[3] * x[, 3])
	J <- matrix(1, 15, 3)
	J[, 2] <- .Uminus(x[, 1]) * x[, 2] * (t[2] * x[, 2] + t[3] * x[, 3])^-2
	J[, 3] <- .Uminus(x[, 1]) * x[, 3] * (t[2] * x[, 2] + t[3] * x[, 3])^-2
	return(f, J)
}
"beale"<-
function(b)
{
	f <- rep(0, 3)
	f[1] <- 1.5 - b[1] * (1 - b[2])
	f[2] <- 2.25 - b[1] * (1 - b[2]^2)
	f[3] <- 2.625 - b[1] * (1 - b[2]^3)
	J <- matrix(c((1 - b[2]), (1 - b[2]^2), (1 - b[2]^3), .Uminus(b[1]), -2 *
		b[1] * b[2], -3 * b[1] * b[2]^2), 3, 2)
	return(f, J)
}
"biggs"<-
function(b)
{
	t <- (1:13) * 0.10000000000000001
	y <- exp(.Uminus(t)) - 5 * exp(-10 * t) + 3 * exp(-4 * t)
	f <- y - b[3] * exp(.Uminus(t) * b[1]) + b[4] * exp(.Uminus(t) * b[2]) - 
		b[6] * exp(.Uminus(t) * b[5])
	J <- matrix(0, 13, 6)
	J[, 1] <- .Uminus(t) * b[3] * exp(.Uminus(b[1]) * t)
	J[, 2] <- t * b[4] * exp(.Uminus(b[2]) * t)
	J[, 3] <- exp(.Uminus(b[1]) * t)
	J[, 4] <- .Uminus(exp(.Uminus(b[2]) * t))
	J[, 5] <- .Uminus(t) * b[6] * exp(.Uminus(b[5]) * t)
	J[, 6] <- exp(.Uminus(b[5]) * t)
	return(f, J)
}
"brown"<-
function(b)
{
	t <- (1:20)/5
	f <- (b[1] + t * b[2] - exp(t))^2 + (b[3] + b[4] * sin(t) - cos(t))^2
	J <- matrix(0, 20, 4)
	J[, 1] <- -2 * (b[1] + t * b[2] - exp(t))
	J[, 2] <- -2 * t * (b[1] + t * b[2] - exp(t))
	J[, 3] <- -2 * (b[3] + b[4] * sin(t) - cos(t))
	J[, 4] <- -2 * sin(t) * (b[3] + b[4] * sin(t) - cos(t))
	return(f, J)
}
"el1"<-
function(b)
{
	f <- rep(0, 3)
	f[1] <- 10 - (b[1]^2 + b[2])
	f[2] <- 7 - (b[1] + b[2]^2)
	f[3] <- 1 - (b[1]^2 - b[2]^3)
	J <- matrix(c(2 * b[1], 1, 2 * b[1], b[2], 2 * b[2], -3 * b[2]^2,  ), 3,
		2)
	return(f, J)
}
"el2"<-
function(b)
{
	f <- rep(0, 6)
	f[1] <- 1 - (b[1]^2 + b[2]^2 + b[3]^2)
	f[2] <- -4 - (b[1]^2 + b[2]^2 + b[3]^2 - 4 * b[3])
	f[3] <- 1 - (b[1] + b[2] + b[3])
	f[4] <- -1 - (b[1] + b[2] - b[3])
	f[5] <- .Uminus((2 * b[1]^3 + 6 * b[2]^2 + 2 * (5 * b[3] - b[1] + 1)^2)
		)
	f[6] <- .Uminus((b[1]^2 - 9 * b[3]))
	J <- matrix(c(2 * b[1], 2 * b[1], 1, 1, (6 * b[1]^2 + 4 * b[1] - 20 * b[
		3] - 4), 2 * b[1], 2 * b[2], 2 * b[2], 1, 1, 12 * b[2], 0, 2 * 
		b[3], (2 * b[3] - 4), 1, -1, (-20 * b[1] + 100 * b[3] + 20), -9
		), 6, 3)
	return(f, J)
}
"madsen"<-
function(t)
{
	f <- rep(0, 3)
	f[1] <-  - ((t[1]^2 + t[2]^2 + t[1] * t[2]))
	f[2] <-  - (sin(t[1]))
	f[3] <-  - (cos(t[2]))
	J <- matrix(0, 3, 2)
	J[1, 1] <- (2 * t[1] + t[2])
	J[1, 2] <- (2 * t[2] + t[1])
	J[2, 1] <- cos(t[1])
	J[3, 2] <-  - (sin(t[2]))
	return(f, J)
}
"mek.1"<-
function(x, y, kmax = 1000, w, int = T, eps = 9.9999999999999995e-07, beta = 
	0.96999999999999997)
{
	if(int == T)
		x <- cbind(1, x)
	sr <- 10^15
	k <- 1
	while(k <= kmax & sr - crossprod(y, w) > eps) {
		d <- pmin(1 - w, 1 + w)
		z <- lsfit(x, y, d^2, int = F)
		sr <- sum(abs(z$resid))
		k <- k + 1
		s <- z$resid * d^2
		alpha <- max(eps, pmax(s/(1 - w), .Uminus(s)/(1 + w)))
		w <- w + (beta/alpha) * s
	}
	coef <- z$coef
	return(coef, w)
}
"mek.rq"<-
function(x, y, kmax = 1000, w, theta = 0.5, int = T, big = 1e+20, eps = 
	9.9999999999999995e-07, beta = 0.96999999999999997)
{
	if(int == T)
		x <- cbind(1, x)
	yw <- big
	k <- 1
	while(k <= kmax & yw - crossprod(y, w) > eps) {
		d <- pmin(theta - w, 1 - theta + w)
		z <- lsfit(x, y, d^2, int = F)
		yw <- sum(rho.rq(z$resid, theta))
		k <- k + 1
		s <- z$resid * d^2
		alpha <- max(eps, pmax(s/(theta - w),  - s/(1 - theta + w)))
		w <- w + (beta/alpha) * s
	}
	coef <- z$coef
	return(coef, w)
}
"model.step"<-
function(lambda, t0, step, model)
{
	sum(abs(model(t0 + lambda * step)$f))
}
"model.step.rq"<-
function(lambda, t0, step, model, theta)
{
	sum(rho.rq(model(t0 + lambda * step)$f, theta))
}
"nll1"<-
function(model, t, k = 2, eps = 9.9999999999999995e-08, beta = 
	0.96999999999999997)
{
#function to compute nonlinear l1 estimate
#	t is the initial value of the unknown parameter
#model is a user-provided function which returns components
#	f=(f_i (x_i , t)
#	J=(grad f_i )
#
	m <- model(t)
	n <- length(m$f)
	zero <- rep(0, n)
	w <- zero
	snew <- sum(abs(m$f))
	sold <- 1e+20
	while(sold - snew > 9.9999999999999995e-08) {
		z <- mek.1(m$J, m$f, k, w, int = F, eps = eps, beta = beta)
		step <- z$coef
		t0 <- t
		l <- nlminb(start = 1, objective = model.step, lower = 0, upper
			 = 1, model = model, t0 = t, step = step)$parameters
		browser()
		t <- t0 + l * step
		m <- model(t)
		sold <- snew
		snew <- sum(abs(m$f))
		w <- lsfit(m$J, z$w, int = F)$resid
		if(max(abs(w)) >= 1)
			w <- w/(max(abs(w)) + eps)
		print(c(t, l, sum(abs(m$f))))
	}
	t
}
"nlrq"<-
function(model, t, k = 2, theta = 0.5, big = 1e+20, nit.max = 100, eps = 
	9.9999999999999995e-08, beta = 0.96999999999999997)
{
#function to compute nonlinear rq estimate
#	t is the initial value of the unknown parameter
#	model is a user-provided function which returns components
#		f=(f_i (x_i , t)
#		J=(grad f_i )
#	theta is the desired quantile
#	k is the number of Meketon steps per iteration
#	eps and eta are tolerance parameters
#
#function returns
#	coef is the value of the parameter at the solution
#	obj is the value of the objective function at the solution
#	nit is the number of "Meketon steps" taken
	m <- model(t)
	n <- length(m$f)
	w <- rep(0, n)
	snew <- sum(rho.rq(m$f, theta))
	sold <- big
	nit <- 0
	while(sold - snew > eps & nit < nit.max) {
		z <- mek.rq(m$J, m$f, k, w, theta = theta, int = F, eps = eps, 
			beta = beta)
		step <- z$coef
		t0 <- t
		l <- nlminb(start = 1, objective = model.step.rq, lower = 0, 
			upper = 1, model = model, t0 = t, theta = theta, step
			 = step)$parameters
		t <- t0 + l * step
		m <- model(t)
		sold <- snew
		snew <- sum(rho.rq(m$f, theta))
		w <- lsfit(m$J, z$w, int = F)$resid
		w1 <- max(pmax(w, 0))
		if(w1 > theta)
			w <- (w * theta)/(w1 + eps)
		w0 <- max(pmax( - w, 0))
		if(w0 > 1 - theta)
			w <- (w * (1 - theta))/(w0 + eps)
		print(c(t, l, sum(rho.rq(m$f, theta))))
		nit <- nit + 1
	}
	return(coef = t, obj = snew, nit)
}
"osb1"<-
function(b)
{
	y <- c(0.84399999999999997, 0.90800000000000003, 0.93200000000000005, 
		0.93600000000000005, 0.92500000000000004, 0.90800000000000003, 
		0.88100000000000001, 0.84999999999999998, 0.81799999999999995, 
		0.78400000000000003, 0.751, 0.71799999999999997, 
		0.68500000000000005, 0.65800000000000003, 0.628, 
		0.60299999999999998, 0.57999999999999996, 0.55800000000000005, 
		0.53800000000000003, 0.52200000000000002, 0.50600000000000001, 
		0.48999999999999999, 0.47799999999999998, 0.46700000000000003, 
		0.45700000000000002, 0.44800000000000001, 0.438, 
		0.43099999999999999, 0.42399999999999999, 0.41999999999999998, 
		0.41399999999999998, 0.41099999999999998, 0.40600000000000003)
	t <- 10 * ((1:33) - 1)
	f <- y - (b[1] + b[2] * exp(.Uminus(t) * b[4]) + b[3] * exp(.Uminus(t) * 
		b[5]))
	J <- matrix(0, 33, 5)
	J[, 1] <- 1
	J[, 2] <- exp(.Uminus(t) * b[4])
	J[, 3] <- exp(.Uminus(t) * b[5])
	J[, 4] <- .Uminus(t) * b[2] * exp(.Uminus(t) * b[4])
	J[, 5] <- .Uminus(t) * b[3] * exp(.Uminus(t) * b[5])
	return(f, J)
}
"osb2"<-
function(b)
{
	y <- y.osb2
	t <- ((1:65) - 1)/10
	f <- y - (b[1] * exp(.Uminus(t) * b[5]) + b[2] * exp(.Uminus((t - b[9])^
		2) * b[6]) + b[3] * exp(.Uminus((t - b[10])^2) * b[7]) + b[4] * 
		exp(.Uminus((t - b[11])^2) * b[8]))
	J <- matrix(0, 65, 11)
	J[, 1] <- exp(.Uminus(t) * b[5])
	J[, 2] <- exp(.Uminus((t - b[9])^2) * b[6])
	J[, 3] <- exp(.Uminus((t - b[10])^2) * b[7])
	J[, 4] <- exp(.Uminus((t - b[11])^2) * b[8])
	J[, 5] <- .Uminus(t) * b[1] * exp(.Uminus(t) * b[5])
	J[, 6] <- .Uminus(b[2]) * (t - b[9])^2 * exp(.Uminus((t - b[9])^2) * b[
		6])
	J[, 7] <- .Uminus(b[3]) * (t - b[10])^2 * exp(.Uminus((t - b[10])^2) * 
		b[7])
	J[, 8] <- .Uminus(b[4]) * (t - b[11])^2 * exp(.Uminus((t - b[11])^2) * 
		b[8])
	J[, 9] <- 2 * b[2] * b[6] * (t - b[9]) * exp(.Uminus((t - b[9])^2) * b[
		6])
	J[, 10] <- 2 * b[3] * b[7] * (t - b[10]) * exp(.Uminus((t - b[10])^2) * 
		b[7])
	J[, 11] <- 2 * b[4] * b[8] * (t - b[11]) * exp(.Uminus((t - b[11])^2) * 
		b[8])
	return(f, J)
}
"powell"<-
function(b)
{
	f <- rep(0, 4)
	f[1] <- .Uminus((b[1] + 10 * b[2]))
	f[2] <- .Uminus(sqrt(5)) * (b[3] - b[4])
	f[3] <- .Uminus((b[2] - 2 * b[3])^2)
	f[4] <- .Uminus(sqrt(10)) * (b[1] - b[4])^2
	J <- matrix(c(1, 0, 0, 2 * sqrt(10) * (b[1] - b[4]), 10, 0, 2 * (b[2] - 
		2 * b[3]), 0, 0, sqrt(5), -4 * (b[2] - 2 * b[3]), 0, 0, .Uminus(
		sqrt(5)), 0, -2 * sqrt(10) * (b[1] - b[4])), 4, 4)
	return(f, J)
}
"rho.rq"<-
function(u, theta)
{
	theta * pmax(u, 0) + (theta - 1) * pmin(u, 0)
}
"rosen"<-
function(b)
{
	f <- rep(0, 2)
	f[1] <- 10 * (b[2] - b[1]^2)
	f[2] <- 1 - b[1]
	J <- matrix(c(20 * b[1], 1, -10, 0), 2, 2)
	return(f, J)
}
"stack"<-
function(b)
{
	y <- stack.loss
	x <- cbind(1, stack.x)
	f <- y - x %*% b
	J <- x
	return(f, J)
}
"watson"<-
function(b)
{
	u <- (1:29)/29
	f <- rep(0, 31)
	f[1:29] <- 1 - (b[2] + 2 * b[3] * u + 3 * b[4] * u^2) + (b[1] + b[2] * 
		u + b[3] * u^2 + b[4] * u^3)^2
	f[30] <-  - b[1]
	f[31] <- 1 - b[2] + b[1]^2
	J <- matrix(0, 31, 4)
	J[1:29, 1] <-  - ((b[1] + u * b[2] - u^2 * b[3] + u^3 * b[4]))
	J[1:29, 2] <- 1 - u * (b[1] + u * b[2] - u^2 * b[3] + u^3 * b[4])
	J[1:29, 3] <- 2 * u - u^2 * (b[1] + u * b[2] - u^2 * b[3] + u^3 * b[4])
	J[1:29, 4] <- 3 * u^2 - u^3 * (b[1] + u * b[2] - u^2 * b[3] + u^3 * b[4
		])
	J[30,  ] <- c(1, 0, 0, 0)
	J[31,  ] <- c(-2 * b[1], 1, 0, 0)
	return(f, J)
}
"wom"<-
function(t)
{
#Womersley's example -- data in xw for x and yw for y
	x <- xw
	y <- yw
	f <- rep(0, 40)
	f <- log10(y) - pmin(log10(x[, 1]), (t[1] + t[2] * x[, 2]))
	J <- matrix(0, 40, 2)
	nc <- log10(x[, 1]) >= t[1] + t[2] * x[, 2]
	J[nc, 1] <- 1
	J[nc, 2] <- x[nc, 2]
	return(f, J)
}
"wood"<-
function(b)
{
	f <- rep(0, 6)
	f[1] <- -10 * (b[2] - b[1]^2)
	f[2] <- 1 - b[1]
	f[3] <- .Uminus(sqrt(90)) * (b[4] - b[3]^2)
	f[4] <- 1 - b[3]
	f[5] <- 2 * sqrt(10) - sqrt(10) * (b[2] + b[4])
	f[6] <- .Uminus((b[2] - b[4]))/sqrt(10)
	J <- matrix(c(-20 * b[1], 1, 0, 0, 0, 0, 10, 0, 0, 0, sqrt(10), 1/sqrt(
		10), 0, 0, -2 * sqrt(90) * b[3], 1, 0, 0, 0, 0, sqrt(90), 0, 
		sqrt(10), -1/sqrt(10)), 6, 4)
	return(f, J)
}
"x.bard"<-
structure(.Data = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 15, 14, 
	13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8, 7, 6,
	5, 4, 3, 2, 1), .Dim = c(15, 3))
"xw"<-
structure(.Data = c(8064, 8064, 8064, 8064, 8064, 8064, 8064, 8064, 8064, 8064, 
	5448, 5448, 5448, 5448, 5448, 5448, 5448, 5448, 5448, 5448, 1680, 1680, 
	1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 528, 528, 528, 528, 528,
	528, 528, 528, 528, 528, 2.3629489603024574, 2.3629489603024574, 
	2.3629489603024574, 2.3629489603024574, 2.3629489603024574, 
	2.3629489603024574, 2.3629489603024574, 2.3629489603024574, 
	2.3629489603024574, 2.3629489603024574, 2.256317689530686, 
	2.256317689530686, 2.256317689530686, 2.256317689530686, 
	2.256317689530686, 2.256317689530686, 2.256317689530686, 
	2.256317689530686, 2.256317689530686, 2.256317689530686, 
	2.1588946459412783, 2.1588946459412783, 2.1588946459412783, 
	2.1588946459412783, 2.1588946459412783, 2.1588946459412783, 
	2.1588946459412783, 2.1588946459412783, 2.1588946459412783, 
	2.1588946459412783, 2.1132713440405748, 2.1132713440405748, 
	2.1132713440405748, 2.1132713440405748, 2.1132713440405748, 
	2.1132713440405748, 2.1132713440405748, 2.1132713440405748, 
	2.1132713440405748, 2.1132713440405748), .Dim = c(40, 2))
"y.bard"<-
c(0.14000000000000001, 0.17999999999999999, 0.22, 0.25, 0.28999999999999998, 
	0.32000000000000001, 0.34999999999999998, 0.39000000000000001, 0.37, 
	0.57999999999999996, 0.72999999999999998, 0.95999999999999996, 
	1.3400000000000001, 2.1000000000000001, 4.3899999999999997)
"y.osb2"<-
c(1.3660000000000001, 1.1910000000000001, 1.1120000000000001, 
	1.0129999999999999, 0.99099999999999999, 0.88500000000000001, 
	0.83099999999999996, 0.84699999999999998, 0.78600000000000003, 
	0.72499999999999998, 0.746, 0.67900000000000005, 0.60799999999999998, 
	0.65500000000000003, 0.61599999999999999, 0.60599999999999998, 
	0.60199999999999998, 0.626, 0.65100000000000002, 0.72399999999999998, 
	0.64900000000000002, 0.64900000000000002, 0.69399999999999995, 
	0.64400000000000002, 0.624, 0.66100000000000003, 0.61199999999999999, 
	0.55800000000000005, 0.53300000000000003, 0.495, 0.5, 
	0.42299999999999999, 0.39500000000000002, 0.375, 0.372, 
	0.39100000000000001, 0.39600000000000002, 0.40500000000000003, 
	0.42799999999999999, 0.42899999999999999, 0.52300000000000002, 
	0.56200000000000006, 0.60699999999999998, 0.65300000000000002, 
	0.67200000000000004, 0.70799999999999996, 0.63300000000000001, 
	0.66800000000000004, 0.64500000000000002, 0.63200000000000001, 
	0.59099999999999997, 0.55900000000000005, 0.59699999999999998, 0.625, 
	0.73899999999999999, 0.70999999999999996, 0.72899999999999998, 
	0.71999999999999997, 0.63600000000000001, 0.58099999999999996, 
	0.42799999999999999, 0.29199999999999998, 0.16200000000000001, 
	0.098000000000000004, 0.053999999999999999)
"yw"<-
c(8064, 8064, 8064, 8064, 8064, 8064, 8064, 8064, 8064, 8064, 5448, 5448, 5448, 
	5196, 4860, 3780, 3542, 3444, 2772, 1764, 1680, 1680, 1680, 1680, 1680, 
	1440, 1344, 1344, 408, 408, 528, 528, 528, 528, 528, 504, 504, 504, 408,
	408)
