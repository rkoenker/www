# remember to activate a graphic device driver
if(is.na(dev.list()*1))
	stop("Please activate a graphic device and try again.")
qss(schuette.x,schuette.y,lambda=60)->qss.o
z_seq(min(schuette.x),max(schuette.y),,100)
qss.fit(schuette.x,qss.o$coef,z)->qss.fit.o
plot(schuette.x,schuette.y)
lines(z,qss.fit.o)
