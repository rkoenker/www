load("Data.d")
attach(D)
s <- abs(y) < 1  # try plotting only the subsample as in RD' s plots
xx <- x[s]
yy <- y[s]
pdf("fig1.pdf",height = 4, width = 7)
par(mfrow=c(1,2),cex = .6)
plot(x,y,xlab = "x", ylab = "y")
plot(xx,yy,xlab = "x", ylab = "y")
dev.off()
