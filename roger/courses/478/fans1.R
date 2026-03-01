#Reproduce Fig 2.1 from Therneau and Grambsch
library(survival5)
d_read.table("fans.dat",header=T)
fans_d[,1]
cens_d[,2]
ps.options(width=6,height=4,horizontal=F)
#ps.options(paper="special",width=6,height=4,horizontal=F)
postscript(file="fans1.ps")
plot(survfit(Surv(fans,cens),type="fleming-harrington"),fun="cumhaz")
#exponential fit
lambda_sum(cens)/sum(fans)
abline(coef=c(0,lambda),lty=2)
dev.off()
