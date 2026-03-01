#Analysis of the AR(1) Melbourne Temperature Example
library(quantreg)
library(splines)
source("melbourne.R")
fitit <- TRUE
plotit <- TRUE
x <- mel.temp[-3650]
y <- mel.temp[-1]
s <- (x<40)
x <- x[s]
y <- y[s]
postscript("fig3aR.ps",horizontal=FALSE,width=6.5,height=6.5)
par(cex=1,pty="s")
xlab <- "yesterday's max temperature"
ylab <- "today's max temperature"
plot(x,y,pch=".",xlab=xlab,ylab=ylab)
taus <- 1:99/100
X <- bs(x,knots=quantile(x),intercept=TRUE)
X[is.na(X)] <- 0
z <- seq(min(x)+5,max(x)-5,length=200)
Z <- bs(z,knots=quantile(x),intercept=TRUE)
fit <- matrix(0,ncol(X),length(taus))
k <- 0
if(fitit){
	for(tau in taus){
		k <- k + 1
		fit[,k] <- rq(y~X-1,tau,ci=FALSE)$coef
		lines(z,Z%*%fit[,k])
		}
	}
if(plotit){
k.plot <- (1:19)*5
	for(k in k.plot){
		lines(z[20:180],Z[20:180,]%*%fit[,k])
		}
}

abline(c(0,1),lty=3)
dev.off()
if(TRUE){
postscript("fig4R.ps",width=7.5,height=5,horizontal=FALSE)
par(cex=1,pty="s")
par(mfrow=c(2,3))
xlab <- "today's max temperature"
for(i in 1:6){
	k <- (i/7)*length(z)
	qtl <- Z[k,]%*%fit
	qtl <- qtl[-1]
	fhat <- akj(qtl,qtl,diff(taus),h=1)$dens
	plot(qtl,fhat,type="l",xlab=xlab,ylab="density")
	title(paste("Yesterday's Temp", format(round(z[k]))))
	}
dev.off()
}
