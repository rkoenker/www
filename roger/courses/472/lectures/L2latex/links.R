#plots of link functions for box-cox models
pdf("links.pdf",width=6.5,height=6.5)
lams <- c(-1,.5,2,4)
eps <- .2
u <- (1:500)/100
plot(u,u-1,type="l",axes=T, xlab="",ylab="")
text(5,4-eps,paste(1))
lines(u,log(u))
text(5,log(5)-eps,paste(0))

for(lam in lams){
	lines(u,(u^lam -1)/lam,lty=2)
	text(5,(5^lam -1)/lam-eps,paste(lam))
	}
text(3,3.5,paste(2))
text(2,3,paste(4))
dev.off()
