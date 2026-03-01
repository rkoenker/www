#Analysis of mc1
#require(logspline)
#par(mfrow=c(2,3))
#for(i in 1:2){
# for(j in 1:3){
# logspline.plot(logspline.fit(S[,i,j]))
# }}

#convert test statistics to decisions
	R <- 1000
	alphas <- c(1,.95,.9,.85)
	ns <- c(100,200)

D <- array(0,c(length(ns),length(alphas),R,2,3))
D[,,,1,] <- (S[,,,1,]< critval.hansen(S[,,,2,]))
D[,,,2,] <- (S[,,,3,]< qnorm(.05))
Dtab <- apply(D,c(1,2,4,5),mean)
