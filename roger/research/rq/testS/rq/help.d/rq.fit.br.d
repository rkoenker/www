.if \n(.g \{\ 
.warn 0
.fp 5 B
.ftr C Bn\}
.BG
.FN  rq.fit.br 
.TL
Quantile Regression Fitting by Exterior Point Methods
.DN
This function controls the details of QR fitting by the simplex approach
embodied in the algorithm of Koenker and d'Orey based on the median
regression algorithm of Barrodale and Roberts.  Typically, options
controlling the construction of the confidence intervals would be passed
via the ... argument of rq().
.CS
rq.fit.br(x, y, tau=0.5, alpha=0.1, ci=T, iid=T, interp=T, tcrit=T)
.RA
.AG x
the design matrix
.AG y
the response variable
.AG tau
the quantile desired, if tau lies outside (0,1) the whole process is estimated.
.OA
.AG alpha
the nominal coverage probability for the confidence intervals 
.AG ci
logical flag if T then compute confidence intervals for the parameters
using the rank inversion method of Koenker (1994).  See rq() for more
details.  If F then return only the estimated coefficients.  Note that
for large problems the default option ci = T can be rather slow.
Note also that rank inversion only works for p>1, an error message is
printed in the case that ci=T and p=1.
.AG iid
logical flag if T then the rank inversion is based on an assumption of
iid error model, if F then it is based on an nid error assumption.
See Koenker and Machado (1999) for further details on this distinction.
.AG interp
As with typical order statistic type confidence intervals the test
statistic is discrete, so it is reasonable to consider intervals that
interpolate between values of the parameter just below the specified
cutoff and values just above the specified cutoff.  If interp =F then
the 2 "exact" values above and below on which the interpolation would
be based are returned.
.AG tcrit
Logical flag if T -  Student t critical values are used, if F then normal
values are used.
.nr XX 0
.ft 1
.fi
.in .75i
.sp
.ti 0
.ta .75i 1i 1.25i 1.5i 1.75i 2i
VALUE:\t\c
.br
Returns an object of class rq for tau in (0,1), or else of class rq.process.
See rq.object and rq.process.object for further details.
.SE

.DT
If tau lies in (0,1) then an object of class rq is returned with various
related inference apparatus.  If tau lies outside [0,1] then an object
of class rq.process is returned.  In this case parametric programming
methods are used to find all of the solutions to the QR problem for
tau in (0,1), the p-variate resulting process is then returned as the
array sol containing the primal solution and dsol containing the dual
solution.  There are roughly O(nlogn) distinct solutions, so users should
be aware that these arrays may be large and somewhat time consuming to
compute for large problems.
.SH REFERENCES 
Koenker, R. and J.A.F. Machado, (1999) Goodness of Fit and Related Inference 
Processes for Quantile Regression,
\fI J. of Am Stat. Assoc.\fR, forthcoming
.SA
rq, rq.fit.fn
.EX
rq.fit.br(stack.x,stack.loss,tau=.73,interp=F)
.KW
quantile regression, simplex method, exterior point, rank test,
.SH Documentation Class
function
.WR
