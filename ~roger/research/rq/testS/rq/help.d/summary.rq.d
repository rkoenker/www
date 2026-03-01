.if \n(.g \{\ 
.warn 0
.fp 5 B
.ftr C Bn\}
.BG
.FN  summary.rq 
.TL
Summary method for Quantile Regression
.DN
Returns a summary list for a quantile regression fit.  A null value
will be returned if printing is invoked.
.CS
summary.rq(object, se="nid", covariance=T)
.RA
.AG object
This is an object of class rq produced by a call to rq().
.OA
.AG se
specifies the method used to compute standard standard errors.  There
are currently three available methods:  
 

   1.  "iid" which presumes that the errors are iid and computes
               an estimate of the asymptotic covariance matrix as in KB(1978).
 
   2.  "nid" which presumes local (in tau) linearity (in x) of the
               the conditional quantile functions and computes a Huber
               sandwich estimate using a local estimate of the sparsity.
 
   3.  "ker" which uses a kernel estimate of the sandwich as proposed
               by Powell(1990).
 

.AG covariance
logical flag to indicate whether the full covariance matrix of the 
estimated parameters should be returned. 
.nr XX 0
.ft 1
.fi
.in .75i
.sp
.ti 0
.ta .75i 1i 1.25i 1.5i 1.75i 2i
VALUE:\t\c
.br
a list is returned with the following components
.RC coefficients
a p by 4 matrix consisting of the coefficients, their estimated standard
errors, their t-statistics, and their associated p-values.
.RC cov
the estimated covariance matrix for the coefficients in the model,
provided that cov=T in the called sequence.
.RC Hinv
inverse of the estimated Hessian matrix returned if cov=T and se!=iid.
.RC J
Outer product of gradient matrix returned if cov=T and se!=iid.
The Huber sandwich is cov = Hinv%*%J%*%Hinv.
.SE

.DT

.SH REFERENCES 
Koenker, R. (2000) \fIQuantile Regression\fR.
.SA
rq, summary.rq.process
.EX
y <- stack.loss
x <- stack.x
summary(rq(y~x,method="fn")) #Compute se's for fit using "nid" method.
summary(rq(y~x,ci=F),se="ker") #default "br" alg, and compute kernel method se's
.KW
quantile regression, confidence intervals, covariance matrix
.SH Documentation Class
function
.WR
