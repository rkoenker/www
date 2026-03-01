# Test of Bennett Model estimation by logit and rq

rdf <- function(x,y, yc){
   fit <- NULL
   for(i in 1:length(yc)){
      fit <- cbind(fit,glm((y > yc[i]) ~ x,
         family = binomial(link = "logit"))$coef)
      }
   fit
   }
      
require(quantreg)
set.seed(1968)
b <- 1
R <- 1000
taus <- 1:4/5
ns <- c(100,500,1000,5000)
A <- array(0,c(4,4,length(ns),R))
B <- rbind(qlogis(taus),rep(1,4))
for(j in 1:length(ns)){
   n <- ns[j]
   x <- rnorm(n)
   for(i in 1:R){
      y <- x * b + rlogis(n)
      b0 <- rq(y ~ x,tau = taus)$coef - B
      b1 <- rdf(x,y,yc = -qlogis(taus)) - B
      A[,,j,i] <- rbind(b0,b1)
      }
   }
save(A,file = "sim1.Rda")
