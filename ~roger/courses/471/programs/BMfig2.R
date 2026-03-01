# plot Boston Marathon Times
#Read Data
x <- read.table("bm.d")
y <- 60 * x[,1] + x[,2] + x[,3]/60
ym <- ts(y,1897)
tm <- ts(1897+1:length(y),1897)
x <- read.table("bw.d")
y <- 60 * x[,1] + x[,2] + x[,3]/60
y <- y[-(1:6)]
yw <- ts(y,1972)
tw <- ts(1972+1:length(y),1972)

#plot(ym, type="n",ylim=c(120,200))
plot(ym, type="n",xlim = c(1895,2020),ylim=c(120,180))
abline
points(yw,col="blue")
points(ym,col="red")
abline(lm(yw ~ tw),col="blue")
abline(lm(ym ~ tm),col="red")
legend(1900,135,c("men","women"),col=c("red","blue"),pch = rep(1,2),lty=rep(1,2))
title("Boston Marathon Winning Times")
