.if \n(.g \{\ 
.warn 0
.fp 5 B
.ftr C Bn\}
.BG
.FN  plot.table.rq 
.TL
Plot Table of Quantile Regression Results
.DN
The function makes an array of p - plots based on an array produced
by table.rq of quantile regression results.  The plots each represent
one parameter of the model specified in the formula argument to
table.rq; the plots consist of the point estimates betahat(tau)
plotted against the taus specified in the table.rq command with
a confidence band as produced by rq().
.CS
plot.table.rq(object,nrow=3,ncol=2)
.RA
.AG object
object of class table.rq containing the array to be plotted.
.OA
.AG nrow
number of rows per page of plots
.AG ncol
number of columns per page of plots
.nr XX 0
.ft 1
.fi
.in .75i
.sp
.ti 0
.ta .75i 1i 1.25i 1.5i 1.75i 2i
VALUE:\t\c
.br
NULL
.SE
plots an array of figures on the current graphics device.
.DT
See rq() and rq.fit.br() table.rq() for further details on control of data.
Obviously, further plotting parameters could be added in a more full-blown
version.  This version is meant just to be illustrative.
.SA
rq,rq.fit.br,table.rq
.EX

.SH Documentation Class
function
.WR
