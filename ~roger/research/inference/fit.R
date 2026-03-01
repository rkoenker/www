formula <- log(dur)~treatment+female+black+hispanic+ndependents+
		factor(quarter)+ recall+young+old+durable+lusd
K3_rq.test.khmal(formula,data=penn46)
K4_rq.test.khmal(formula,data=penn46,location.scale=F)
postscript("L.ps",horizontal=F)
plot(K4,ncol=4)
dev.off()
postscript("LS.ps",horizontal=F)
plot(K3,ncol=4)
dev.off()
