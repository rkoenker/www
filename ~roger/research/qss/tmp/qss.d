.BG
.FN qss
.TL
Quantile Smoothing Splines 
.DN
Compute the coefficients of the nonparametric quantile smoothing splines 
.CS
qss(x, y, w=rep(1, length(y)), theta=0.5, lambda=-1, method = "bc", penalty = "Loo",
constraint="none", coef=rep(0, nvar), pswitch=F, maxiter=500, penpar=1/lam^0.5,
tmin = 1/(length(y)+nh), kmax = sqrt(length(y)))
.RA
.AG x
vector (matrix) of independent variable(s).
.AG y
vector (matrix) of dependent variable(s).
.OA
.AG w
vector of weights summing to 1, one for each observation, default to 1/length(x).
.AG theta
desired $ theta$th quantile smoothing spline.  0 <= theta <= 1.  If theta < 0
or >1, the whole spectrum of quantile smoothing splines in theta will be 
returned. Default to the median quantile smoothing spline.
.AG lambda
smoothing parameter; $lambda = 0$ corresponses to interpolating splines while
$lambda = \infty$ corresponses to linear regression quantile.  If lambda < 0,
the whole spectrum of solution in lambda will be returned and is the default.
.AG method
"br" uses the algorithm of Barrodale and Roberts (1973) while "bc" uses
the algorithm of Bartels and Conns (1980).  Default to "bc".
.AG penalty
"L1" uses the $L sub 1$ roughness penalty while "Loo" uses the $L sub \infty$ roughness penalty.
.AG constraint
"increase" imposes positive monotonicity, "decrease" imposes negative monotonicity, 
"concave" imposes concavity and "convex" imposes convexity.  Default to "none".
.AG coef
vector of intial coefficients, default to zeros.
.AG pswitch
control printing of intermediate steps.  No printing by default.
.AG maxiter
maximum number of iteration allowed.  500 by default.
.AG penpar
penalty parameter for "bc", default to 1/sqrt(lambda).
.AG tmin
the smallest value of theta at which the parametric programming should 
start; default to 1/(2*length(y)).
.AG kmax
the largest effective dimensionality at which the parametric programming in
lambda should stop; default to sqrt(length(y)).
.RT If theta is in [0,1] and lambda >= 0
.AG coef
coefficients for the quantile smoothing spline.  
See qss.fit for detail.
.AG resid
vector of residuals.
.AG exit 
exit code (method="br"): 0 -- probably nonunique solution; 1 -- OK
.AG ifl 
exit code (method="bc"):  1 -- OK;
2 -- unsuccessful termination, constraints cannot be satisfied, problem is infeasible;
3 -- limit imposed by maxiter reached without finding a solution;
4 -- program aborted, numerical difficulties, due to ill-conditioning.
.AG niter
number of cycles taken.
.AG fidel
fidelity of the objective function.
.AG rough
roughness of the objective function.
.AG obj
value of the objective function; obj = fidel + lambda * rough.
.AG k
effective dimensionality of the fitted quantile.
.AG penalty
penalty type.
.RT If theta is not in [0,1] or lambda < 0
.AG sol
a (nvar+6) by nsol matrix whose first row contains the theta `breakpoints' 
of the quantile smoothing splines, the second row contains the lambda
`breakpoints', the third row the exit code, the fourth row
the fidelity, the fifth row the roughness and the sixth row the effective
dimensionality of the fits corresponding to the theta and lambda in every
column.  The remaining rows contain the coefficients of the quantile
smoothing splines for the range of theta and lambda from that column
to the next.
.SH REFERENCES
Barrodale, I. and F.D.K. Roberts, "An Improved Algorithm for Discrete $l sub 1$
Linear Approximation", SIAM J. Numer. Anal., 10, 1973, 839-848.

Bartels, R. and A. Conn, "Linearly Constrained Discrete $l sub 1$ Problem",  ACM
Transactions on Mathematical Software, 6, 1980, 594-608.
.SA
qss.fit
.EX
.KW ~keyword
qss(corn.rain,corn.yield,lambda=40,penalty="Loo")->qss.o
z_seq(min(corn.rain),max(corn.rain),,100)
qss.fit(corn.rain,qss.o$coef,z)->qss.fit.o
plot(corn.rain,corn.yield)
lines(z,qss.fit.o,type="l")
.WR
