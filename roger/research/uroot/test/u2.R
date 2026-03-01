"uroot.test" <-
function(y, p=1, l=4, trend = FALSE, res = 0, 
	scores = c("wilcoxon", "normal", "sign"))
{
#Compute Hasan-Koenker Regression Rankscore Test of the Unit Root Hypothesis
#y is a univariate timeseries
#p is the lag length of the Augmented Dickey-Fuller Model
#l is the lag truncation length for the Newey-West Estimates
#if trend=T then a trend term is included in the ADF model
#score should be one of the following:
#	"wilcoxon" for the logistic model
#	"normal" for the Gaussian model
#	"sign" for the Laplace model
#res controls the treatment of zero residuals in the estimation of Sigma
#
	M <- make.adf.model(y, p, trend = trend)
	dy <- M$dy
	z <- M$z
	x0 <- M$x0
	fit <- rq(dy~x0, tau = -1)
	x1hat <- as.matrix(qr.resid(qr(cbind(1, x0)), z))
	sgmzz <- (crossprod(x1hat))/((length(dy)^2))
	Stats <- matrix(0,3,length(scores))
	for(i in 1:length(scores)){
	   R <- GJKP.ranks(fit, scores[i])
	   S <- (t(x1hat) %*% R)/sqrt(crossprod(x1hat))
	   adf.coef <- lm(dy ~ cbind(z, x0))$coef[2]
	   w <- uroot.nuisance(x0, dy, l, res, trend, scores[i])
	   delta2 <- w$Sigma[1,2]^2/(w$Sigma[1,1]*w$Sigma[2,2])
	   D <- determinant(w$Sigma)
	   w1 <- sqrt(w$Sigma[1, 1]/D)
	   w2 <- (w$Sigma[1, 2] * sqrt(sgmzz))/(sqrt(D * w$Sigma[1, 1]))
	   w3 <- (w$Sigma[1, 1] - w$sgmu2)/(2 * sgmzz)
	   rtest <- c(w1 * S - w2 * (length(dy) * adf.coef - w3))
	   Stats[,i] <- c(S,delta2,rtest)
	  }
	#rescale the unmodified test statistics
	Stats[1,] <- Stats[1,]/sqrt(c(1/12,1,1/4))
	return(Stats)
}
