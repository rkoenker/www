.if \n(.g \{\ 
.warn 0
.fp 5 B
.ftr C Bn\}
.BG
.FN  rq.fit.fn 
.TL
Quantile Regression Fitting via Interior Point Methods
.DN
This is a lower level routine called by rq() to compute quantile
regression methods using the Frisch-Newton algorithm.
.CS
rq.fit.fn(x, y, tau=0.5, int=F, beta=0.99995, eps=1e-06)
.RA
.AG x
The design matrix
.AG y
The response vector
.OA
.AG tau
The quantile of interest, must lie in (0,1)
.AG int
logical flag, if T an intercept should be appended to x, if not, not.
.AG beta
technical step length parameter -- alter at your own risk!
.AG eps
tolerance parameter for convergence.  In cases of multiple optimal solutions
there may be some descrepancy between solutions produced by method "fn"
and method "br".  This is due to the fact that "fn" tends to converge to
a point near the centroid of the solution set, while "br" stops at a
vertex of the set.  
.nr XX 0
.ft 1
.fi
.in .75i
.sp
.ti 0
.ta .75i 1i 1.25i 1.5i 1.75i 2i
VALUE:\t\c
.br
returns an object of class rq, which can be passed to summary.rq to obtain
standard errors, etc.
.SE
none
.DT
The details of the algorithm are explained in Koenker and Portnoy (1997).
The basic idea can be traced back to the log-barrier methods proposed by
Frisch in the 1950's for constrained optimization.  But the current
implementation is based on proposals by Mehrotra and others in the
recent (explosive) literature on interior point methods for solving linear 
programming problems.  This version of the algorithm is designed for
fairly large problems, for very large problems see rq.fit.pfn.
.SH REFERENCES 
Koenker, R. and S. Portnoy (1997).
The Gaussian Hare and the Laplacian Tortoise:
Computability of Squared-error vs. Absolute-error Estimators, with discussion,
\fI Statistical Science\fR, 12, 279-300.
.SA
rq, rq.fit, rq.fit.br, rq.fit.pfn
.EX

.KW
quantile regression, linear programming, interior point methods
.SH Documentation Class
function
.WR
