load("flies.Rdata")
attach(flies)
density <- begin -sum(begin*num)/sum(num)
pmale <- prbegin -sum(prbegin*num)/sum(num)
sized <- size -sum(size*num)/sum(num)
formula <- log(age)~ factor(female) + sized +density + pmale 
fitA <- summary(lm(formula,weight=num))$coefficients
dnames <- dimnames(fitA)
p <- nrow(fitA)
taus <- c(.01,1:19/20,96:99/100,995:999/1000)
fitA <- array(fitA,c(p,4,length(taus)))
dimnames(fitA) <- list(dnames[[1]],dnames[[2]],paste(taus))
for(i in 1:length(taus)){
	print(taus[i]) 
	f <- rq(formula, taus[i], weight=num, method="fn",eps=1e-5)
	fitA[,,i] <- summary(f)$coefficients
	}
