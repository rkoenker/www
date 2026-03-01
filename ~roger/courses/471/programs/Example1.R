require(tseries)
require(zoo)

ibm <- get.hist.quote("ibm",start="1998-01-01",quote="AdjClose")

fom <-function(x) as.Date(sub("..$","01",format(x)))
ibmm <- aggregate(ibm,fom(time(ibm)),head,1)
spc <- get.hist.quote("^gspc",start="1998-01-01",quote="Close")

spcm <- aggregate(spc,fom(time(spc)),head,1)


plot(coredata(merge(spcm,ibmm)))
f <- lm(ibmm ~ spcm)
abline(f)
summary(f)


