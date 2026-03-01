.if \n(.g \{\ 
.warn 0
.fp 5 B
.ftr C Bn\}
.BG
.FN rq.object
.TL
Linear Quantile Regression Process Object 
.DN
These are objects of class\fBrq.process. \fP
They represent the fit of a linear conditional quantile function model. 
.SH  GENERATION 
This class of objects is returned from the \fBrq\fP
 function 
to represent a fitted linear quantile regression model. 
.SH  METHODS 
The \fB"rq.process"\fP class of objects has 
methods for the following generic 
functions: 
\fBeffects\fP, \fBformula\fP
, \fBlabels\fP
,  \fBmodel.frame\fP
, \fBmodel.matrix\fP
, \fBplot\fP
, \fBpredict\fP
, \fBprint\fP
, \fBprint.summary\fP
, \fBsummary\fP
,  

.SH  STRUCTURE 
The following components must be included in a legitimate \fBrq.process\fP
 object. 
.nr XX 0
.ft 1
.fi
.in .75i
.sp
.ti 0
.ta .75i 1i 1.25i 1.5i 1.75i 2i
VALUE:\t\c
.br
.RC sol
The primal solution array.  This is a (p+3)  by  J  matrix  whose  
first  row  contains  the 'breakpoints'   tau_1,tau_2,...tau_J,   
of   the  quantile function, i.e. the values in [0,1] at which  the  
solution changes,  row  two  contains  the  corresponding quantiles 
evaluated at the mean design point, i.e. the inner product of  
xbar  and  b(tau_i), the third row contains the value of the objective
function evaluated at the corresponding tau_j, and the last p rows 
of the matrix give b(tau_i).  The solution b(tau_i) prevails from  
tau_i to tau_i+1.  Portnoy (1991) shows that J=O_p(n log n).
.RC dsol
The dual solution array.  This is a
n by J matrix  containing the  dual  solution  corresponding to sol,
the ij-th entry is 1 if y_i > x_i b(tau_j), is 0 if y_i <  x_i
b(tau_j),   and  is between 0 and 1 otherwise, i.e. if the
residual is zero.  See  Gutenbrunner  and  Jureckova(1991)
for    a    detailed   discussion   of   the   statistical
interpretation of dsol.  The use of dsol in inference is described
in Gutenbrunner, Jureckova, Koenker, and Portnoy (1994).
.DT
These arrays are computed by parametric linear programming methods
using using the exterior point (simplex-type) methods of the 
Koenker-d'Orey algorithm based on Barrodale and Roberts median
regression algorithm.
.SH  REFERENCES 


[1] Koenker, R.W. and Bassett, G.W. (1978). Regression quantiles,
\fIEconometrica\fR, 46, 33-50.
[2] Koenker, R.W. and d'Orey (1987,1994). Computing Regression Quantiles.
\fIApplied Statistics\fR, 36, 383-393, and 43, 410-414.
[3] Gutenbrunner, C. Jureckova, J. (1991).
Regression quantile and regression rank score process in the
linear model and derived statistics, \fIAnnals of Statistics\fR, 20, 305-330.
[4] Gutenbrunner,  C.,  J.  Jureckova,   Koenker,  R.  and
Portnoy,  S.(1994)  "Tests  of Linear Hypotheses  based on Regression  
Rank   Scores",   \fIJournal   of   Nonparametric Statistics\fR, 
(2), 307-331.
[5]  Portnoy, S. (1991).  Asymptotic behavior of the number of regression
quantile breakpoints, \fI SIAM Journal of  Scientific
and  Statistical Computing\fR, 12, 867-883.
.SA
\fBrq\fP.
.KW
classes
regression
methods
.SH Documentation Class
function
.WR
