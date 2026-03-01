# Make tables for sim2
require(quantreg)
ns <- c(100,500,1000,5000)
taus <- 1:4/5
load("sim2.Rda")
eves <- (1:8) *2
Bias <- matrix(aperm(apply(A,1:3,mean),c(1,3,2)),16,4)
RMSE <- matrix(aperm(sqrt(apply(A^2,1:3,mean)),c(1,3,2)),16,4)
aBias <- Bias[-eves,]
bBias <- Bias[eves,]
aRMSE <- RMSE[-eves,]
bRMSE <- RMSE[eves,]
eves <- (1:4) *2
nord <- c(eves - 1,eves)
scale <- rep(sqrt(ns),2)
Bias <- cbind(aBias[nord,]*scale, bBias[nord,]*scale)
RMSE <- cbind(aRMSE[nord,]*scale, bRMSE[nord,]*scale)
RMSE21 <- RMSE[5:8,1:4]
RMSE12<- RMSE[1:4,5:8]
RMSE[1:4,5:8] <- RMSE21
RMSE[5:8,1:4] <- RMSE12
load("lim2.Rda") # asymptotic se's
tab <- cbind(rbind(RMSE[1:4,1:4],A[3,],RMSE[5:8,1:4],A[4,]),
	rbind(RMSE[1:4,5:8],A[1,],RMSE[5:8,5:8],A[2,]))
dimnames(tab) <- list(
	rep(c(paste("$n =$",ns),paste("$n = \\infty$")),2), 
	rep(paste("$\\tau = $",taus),2))
rgrp <- c("$\\alpha$","$\\beta$")
cgrp <- c("Quantile Regression", "Distributional Regression")
cap <- "Root Mean Squared Error in 1000 Replications: Cox Model"
lab <- "sim2"
latex.table(tab, file = "sim2", rowlabel = "", dec = 3,
        rgroup = rgrp, cgroup = cgrp, caption = cap,label = lab)


