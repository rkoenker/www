# Small timing experiment to compare kuantile and quantile

require(kuantile)
set.seed(1446)

ns <- 10^(2:7)
R <- 10
T <- array(NA,c(R,length(ns),5)) 
eps <- 20 * .Machine$double.eps
for(j in 1:length(ns)){
    for(i in 1:R){
	y <- rnorm(ns[j])
	T[i,j,1] <- system.time(qy <- quantile(y,.5))[1]
	T[i,j,2] <- system.time(ky <- kuantile(y,.5))[1]
	stopifnot(abs(qy - ky) < eps)
	T[i,j,3] <- system.time(qy <- quantile(y))[1]
	T[i,j,4] <- system.time(ky <- kuantile(y))[1]
	stopifnot(abs(qy - ky) < eps)
	T[i,j,5] <- system.time(sort(y))[1]
	}
    }

Tab <- apply(T,2:3,mean)
dimnames(Tab) <- list(paste(ns),c(rep(c("quantile","kuantile"),2),"sort"))
