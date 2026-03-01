.BG
.FN qss.fit
.TL
Quantile Smoothing Splines Fit
.DN
Compute the fitted values of the quantile smoothing splines at the grid z
.CS
qss.fit(x, g, z)
.RA
.AG x
vector (matrix) of independent variable(s).
.AG g
coefficients vector returned by qss
.AG z
vector or list of grid at which fitted values of the quadratic spline are desired;
one component for each independent variable.
.RT
vector (matrix) of fitted values at the grid z
.SA
qss
.EX
qss(evap.x[,c(1,8)],evap.y,lambda=c(50,150),penalty="Loo")->qss.o
z_list(seq(min(evap.x[,1]),max(evap.x[,1]),,25),seq(min(evap.x[,8]),max(evap.x[,8]),,25))
qss.fit(evap.x[,c(1,8)],qss.o$coef,z)->qss.fit.o
persp(qss.fit.o)
.KW ~keyword
.WR
