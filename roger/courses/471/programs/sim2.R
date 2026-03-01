# A Simple Simulation of some least squares confidence ellipses"

#Generate a pair of dependent covariates
X <- matrix(rnorm(200,3,1),100,2)
r <- .7
X <- X %*% matrix(c(1,r,r,1),2,2)

#Setup plot
X11()
plot(1:10,1:10,xlab="x1",ylab="x2",type="n",pty="s")

points(X)
require(ellipse)
lines(ellipse(cov(X),centre=apply(X,2,mean)))

#Now plot some confidence ellipses for regression fits

X11()
plot(1:4,1:4,xlab="beta1",ylab="beta2",type="n",pty="s")
lines(ellipse(solve(cov(X)*100),centre=c(2,2)),col=2)
points(2,2)

for(i in 1:10){
	y <- 2*X[,1] + 2 * X[,2] + rnorm(100)
	f <- lm(y ~ X)
	lines(ellipse(f,which = 2:3),col="grey")
	}
