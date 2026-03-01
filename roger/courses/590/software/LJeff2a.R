# Tables for LJeff.R simulation
require(Hmisc) # For latex
taus <- c(.1, .25,.5,.75,.9)
ns <- c(25,50,100,200)
dfs <- paste("d",1:3,sep = "")
# Three target densities as in Dunson and Taylor
require(nor1mix)
nM2 <- norMix(mu = c(-0.6,-0.5,-0.4,0.5,1),sig2 = c(0.4,0.5,0.6,0.95,1.15)^2)
nM3 <- norMix(mu = c(-0.5,-0.5,-0.5,0.75,1.5),sig2 = c(0.25,0.25,0.25,0.5,1)^2)
qd1 <- function(p) qnorm(p)
qd2 <- function(p) qnorMix(p,nM2)
qd3 <- function(p) qnorMix(p,nM3)
dd1 <- function(x) dnorm(x)
dd2 <- function(x) dnorMix(x,nM2)
dd3 <- function(x) dnorMix(x,nM3)


# Density plot to check D & T's Figure 1.
pdf("denfig.pdf", height = 5, width = 7)
x <- -500:500/100
plot(x,dd1(x),type = "l", ylim = c(0,1), ylab = "f(x)")
lines(x,dd2(x),col = 2)
lines(x,dd3(x),col = 3)
dev.off()


load("LJeff2.Rda")
colnames <- paste("$\\tau = $",round(taus,2))
rownames <- rep(paste("$n = $",ns),2)
rgrp <- c("SampleQ","BayesQ")
Qtrue <- matrix(0,3,length(taus))
Qtrue[1,] <- qd1(taus)
Qtrue[2,] <- qd2(taus)
Qtrue[3,] <- qd3(taus)
rlabel <- paste("Distribution", 1:3)
for(i in 1:3){
   A <- FEJ[[i]]$A
   MSE <- aperm(sqrt(apply((A - Qtrue[i,])^2,1:3,mean)),c(2,3,1))
   MSE <- rbind(MSE[1,,],MSE[2,,])
   dimnames(MSE) <- list(rownames,colnames)
   cap <- paste("Root Mean Square Errors:  Normal Mixture", i, 
       ", 200 replications, 5000 Metropolis-Hastings steps with 
       Jeffreys substitute Likelihood")
   fname <- paste("Jeff", i, ".tex", sep = "")
   latex(MSE,file = fname,rowlabel = rlabel[i], rgroup = rgrp,
	   caption = cap, dec = 3, caption.loc = "bottom")
}
