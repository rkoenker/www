# Test of Cox Model estimation by cloglog and rq

rdf <- function(x,y, yc){
   fit <- NULL
   for(i in 1:length(yc)){
      fit <- cbind(fit,glm((y > yc[i]) ~ x,
         family = binomial(link = "cloglog"))$coef)
      }
   fit
   }
      
require(quantreg)
require(evd)
set.seed(1968)
ns <- c(100,500,1000,5000)
b <- 1
R <- 1000
taus <- 1:4/5
A <- array(0,c(4,4,R))
A <- array(0,c(4,4,length(ns),R))
B0 <- rbind(qgumbel(taus),rep(1,4))
B1 <- rbind(-qgumbel(taus),rep(1,4))
for(j in 1:length(ns)){
   n <- ns[j]
   x <- rnorm(n)
   for(i in 1:R){
      y <- x * b + rgumbel(n)
      b0 <- rq(y ~ x,tau = taus)$coef - B0
      b1 <- rdf(x,y,yc = qgumbel(taus)) - B1
      A[,,j,i] <- rbind(b0,b1)
      }
   }
save(A,file = "sim2.Rda")
