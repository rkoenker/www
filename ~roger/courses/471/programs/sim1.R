# A Simple Simulation of some least squares confidence ellipses"

#Generate a pair of independent covariates
X <- matrix(rnorm(200,5,1),100,2)

#Setup plot
plot(1:10,1:10,xlab="x1",ylab="x2",type="n",pty="s")

points(X)
require(ellipse)
lines(ellipse(cov(X),centre=apply(X,2,mean)))

#Now plot some confidence ellipses for regression fits

plot(1:3,1:3,xlab="beta1",ylab="beta2",type="n",pty="s")
lines(ellipse(solve(cov(X)*100),centre=c(2,2)),col=2)
points(2,2)

for(i in 1:10){
	y <- 2*X[,1] + 2 * X[,2] + rnorm(100)
	f <- lm(y ~ X)
	lines(ellipse(f,which = 2:3),col="grey")
	}
