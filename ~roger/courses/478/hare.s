#fit meyer data for 472 problem with KST-hare-heft methods
library(hare)
library(heft)
#attach(data) #this uses full 3376 observations not the reduced 500 sample in d
attach(d) #this uses reduced sample of 500
x_cbind(ben,wage)
dimnames(x)_list(NULL,c("ben","wage"))

#optional initial fit without covariates to transform timescale
#z_heft.fit(dur,cens)
#tdur_-log(1-pheft(dur,z))
#z_hare.fit(tdur,cens,x)


par(mfrow=c(2,2))
z.hare_hare.fit(dur,cens,x)
hare.summary(z.hare)
xbar_apply(x,2,"mean")
hare.plot(z.hare,cov=xbar,which=0,what="h",time=10)
hare.plot(z.hare,cov=xbar,which=1,what="h",time=10)
x_cbind(ben+rnorm(length(ben))/100,wage)
z.hare_hare.fit(dur+rnorm(len(dur))/10,cens,x)
hare.summary(z.hare)
xbar_apply(x,2,"mean")
hare.plot(z.hare,cov=xbar,which=0,what="h",time=10)
hare.plot(z.hare,cov=xbar,which=1,what="h",time=10)
