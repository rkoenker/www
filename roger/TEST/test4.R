require(quantreg)
test.crq <- read.csv("test.crq.csv",header = TRUE)
test.crq$site_no <- factor(test.crq$site_no)
my<-1990

fit <- crq(Surv(log10(cl),cl_bdl,type="left") ~ I(year.day-my) + log10(q) + sin(2*pi*(I(year.day-my))) + cos(2*pi*(I(year.day-my)))+ site_no + I(year.day-my):site_no + log10(q):site_no + sin(2*pi* (I(year.day-my))):site_no + cos(2*pi*(I(year.day-my))):site_no,data=test.crq,method="PengHuang",contrasts=list(site_no="contr.sum"))

summary.fit <-summary(fit,taus=c(0.10,0.25,0.50,0.75,0.90), se="boot",R=1000)

