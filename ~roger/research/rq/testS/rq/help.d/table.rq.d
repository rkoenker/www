.if \n(.g \{\ 
.warn 0
.fp 5 B
.ftr C Bn\}
.BG
.FN  table.rq 
.TL
Table of Quantile Regression Results
.DN
Function to produce a table of quantile regression results for a group
of specified quantiles.
.CS
table.rq(formula, taus=c(0.05, 0.25, 0.5, 0.75, 0.95), method="br", ...)
.RA

.OA
.AG formula
formula for the linear model, see rq()
.AG taus
quantiles of interest
.AG method
algorithmic method, for large problems method="fn" would be preferred.
.AG ...
other optional arguments passed to rq()
.nr XX 0
.ft 1
.fi
.in .75i
.sp
.ti 0
.ta .75i 1i 1.25i 1.5i 1.75i 2i
VALUE:\t\c
.br
The function returns an array with dimensions (p,m,3), where p is the dimension
of the parameter vector of the model specified by formula, m is the number
of quantiles specified by tau.  For each coefficient at each tau there
is a point estimate and lower and upper limits of a confidence interval.
The object returned is of class table.rq and can be plotted, or formatted
into a latex style table.
.DT
This is only implemented for method="br", but modifications for "fn"
would be quite straightforward.
.SA
rq, rq.fit.br, plot.table.rq, latex.table.rq
.EX
plot(table.rq(stack.loss~stack.x)#plot results of a quantile regression
latex.table.rq(table.rq(stack.loss~stack.x))#make latex table 
.SH Documentation Class
function
.WR
