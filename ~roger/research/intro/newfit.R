# Attempt to reproduce birthweight fit  May 2011
require(quantreg)
load("natality.Rda")
d <- cbind(d, m.wtgain2 = d$m.wtgain^2)
formula <- weight ~ black+married+boy+tri2+tri3+novisit+ed.hs+ed.smcol+
        ed.col+mom.age+smoke+cigsper+m.wtgain+mom.age2+m.wtgain2
fit3.ols <- summary(lm(formula,data = d))$coefficients
p <- nrow(fit3.ols)
taus <- c(1:4/100, 1:19/20)
fit3 <- array(fit3.ols,c(p,4,length(taus)))
for(i in 1:length(taus)){
        print(taus[i])
        f <- rq(formula, taus[i], data = d, method="fn")
        fit3[,,i] <- summary(f)$coefficients
        }


