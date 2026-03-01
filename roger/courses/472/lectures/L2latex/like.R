#this is a source file for the log-likelihood plot for 472 lecture 2.
library("MASS")
pdf("like2.pdf",width=6.5,height=6.5)
#generate the data:
x <- exp(rnorm(50))
y <- exp(1+.5*log(x) +rnorm(50)/3)
z <- boxcox(y~x,plotit=T)
dev.off()
