.if \n(.g \{\ 
.warn 0
.fp 5 B
.ftr C Bn\}
.BG
.FN rq.object
.TL
Linear Quantile Regression Object 
.DN
These are objects of class\fBrq. \fP
They represent the fit of a linear conditional quantile function model. 
.SH  GENERATION 

This class of objects is returned from the \fBrq\fP
 function 
to represent a fitted linear quantile regression model. 
.SH  METHODS 
The \fB"rq"\fP class of objects has methods for the following generic 
functions: 
\fBcoef\fP, \fBeffects\fP
, \fBformula\fP
, \fBlabels\fP
,  \fBmodel.frame\fP
, \fBmodel.matrix\fP
, \fBplot\fP
, \fBpredict\fP
, \fBprint\fP
, \fBprint.summary\fP
, \fBresiduals\fP
, \fBsummary\fP
,  

.SH  STRUCTURE 
The following components must be included in a legitimate \fBrq\fP
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
.RC coefficients
the coefficients of the quantile regression fit. 
The names of the coefficients are the names of the 
single-degree-of-freedom effects (the columns of the 
model matrix). 
If the model was fitted by method "br" with ci=T, then
the coefficient component consists of a matrix whose
first column consists of the vector of estimated coefficients
and the second and third columns are the lower and upper
limits of a confidence interval for the respective coefficients.
.RC residuals
the residuals from the fit. 
.RC contrasts
a list containing sufficient information to construct the contrasts 
used to fit any factors occurring in the model. 
The list contains entries that are either matrices or character vectors. 
When a factor is coded by contrasts, the corresponding contrast matrix 
is stored in this list. 
Factors that appear only as dummy variables and variables in the model 
that are matrices correspond to character vectors in the list. 
The character vector has the level names for a factor or the column 
labels for a matrix. 
.RC model
optionally the model frame, if \fBmodel=T\fP. 
.RC x
optionally the model matrix, if \fBx=T\fP. 
.RC y
optionally the response, if \fBy=T\fP. 

.DT
The coefficients, residuals, and effects may be extracted 
by the generic functions of the same name, rather than 
by the \fB$\fP operator.   
For pure \fBrq\fP objects this is less critical 
than for some of the inheritor classes.  
Note that the extractor function \fBcoef\fP returns a vector with missing values 
omitted.  
.SA
\fBrq\fP,  \fBcoefficients\fP.   
.KW
classes
regression
methods
.SH Documentation Class
function
.WR
