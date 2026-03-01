#This is a new version of the Chicago Land value model using the nprq package
require(nprq)
attach(chiland)
x <- xqtr + .05* runif(1194)
y <- yqtr + .05* runif(1194)
z <- .2 * price/sf
fit <- rqss(z ~ qss(cbind(x,y),lambda = .35, ndum = 10000))
rgl.clear()
plot(fit,rgl=TRUE,ncol=500)
