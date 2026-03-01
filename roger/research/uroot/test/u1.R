"uroot.nuisance" <-
function(x, y, l, res = 0, trend = T, score = "wilcoxon")
{
# x is the adf X matrix under H_o, without intercept term.  
# first p terms are \delta y_t and (p+1)th term is trend if trend==T
# y is the first difference of the univariate time series
# 
# b is the coefficient vector of l_1 fit of the unit root model under H_0 
# r is the residual vector associated with the l_1 fit.  
# p is the number of autoregression terms in the difference
# l is the lag truncation parameter for estimating variance-covariance matrix
# w are the weights for the Newey-West estimator
# res controls the treatment of zero residuals in computing Fhat
#
#Computation of Uhat
	if(is.matrix(x)) n <- nrow(x)
	else n <- length(x)
	fit <- rq(y~x)
	r <- fit$residuals
	b <- fit$coefficients
	ut <- rep(0, n)
	if(trend == T) {
		p <- ncol(x) - 1
		xb <- b[p + 2] * x[, (p + 1)]
		ut <- filter(r + xb, b[2:(p + 1)], "rec")
	}
	else ut <- y - r
	w2 <- c(1/2, 1 - (1:l)/(1 + l))
	z <- ut
	z <- c(rep(0, l), z)
	U <- z
	for(i in 1:l)
		U <- ts.intersect(U, lag(z,  - i))
	sgmu2 <- sum(ut^2)/n
	Sigma <- matrix(0, 2, 2)
	Sigma[1, 1] <- (2 * crossprod(crossprod(U)[1,  ], w2))/n	
	#Computation of empirical distribution function
	if(res == 0)
		Fhat <- (Rank(r) - 0.5)/n
	else Fhat <- (rank(r) - 0.5)/n
	if(score == "wilcoxon")
		vt <- Fhat - 0.5
	else if(score == "normal")
		vt <- qnorm(Fhat)
	else if(score == "sign") {
		eps <- .Machine$double.eps
		vt <- 0.5 * sign(Fhat - 0.5) * (abs(Fhat - 0.5) > eps)
	}
	else stop("Invalid Score Function")
	z <- vt	#Computation of nuisance parameter estimates
	z <- c(rep(0, l), z)
	V <- z
	for(i in 1:l)
		V <- ts.intersect(V, lag(z,  - i))
	w1 <- c(1, 1 - (1:l)/(1 + l))
	w <- 1 - (1:l)/(1 + l)
	suv <- crossprod(crossprod(U, V)[1,  ], w1)/n
	svu <- crossprod(crossprod(V, U)[1, -1], w)/n
	Sigma[1, 2] <- suv + svu
	Sigma[2, 1] <- suv + svu
	Sigma[2, 2] <- (2 * crossprod(crossprod(V)[1,  ], w2))/n
	return(list(Sigma=Sigma,sgmu2=sgmu2))
}
