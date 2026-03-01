.if \n(.g \{\ 
.warn 0
.fp 5 B
.ftr C Bn\}
.BG
.FN  rrs.test 
.TL
Quantile Regression Rankscore Test
.DN
Function to compute regression rankscore test of a linear hypothesis
based on the dual quantile regression process.  A test of  the 
hypothesis, 
is carried out by estimating the restricted model and constructing
a test based on the dual process under the restricted model.  The
details of the test are described in GJKP(1993).  The test has a
Rao-score, Lagrange-multiplier interpretation since in effect it
is based on the value of the gradient of unrestricted quantile regression
problem evaluated under the null.  This function will eventually be
superseded by a more general anova() method for rq.
.CS
rrs.test(x0, x1, y, v, score="wilcoxon")
.RA
.AG x0
the matrix of maintained regressors, a column of ones is appended automatically.
.AG x1
matrix of covariates under test.
.OA
.AG y
response variable, may be omitted if v is provided.
.AG v
object of class rq.process generated e.g. by rq(y~x0,tau=-1)
.AG score
Score function for test (see rq.ranks())
.nr XX 0
.ft 1
.fi
.in .75i
.sp
.ti 0
.ta .75i 1i 1.25i 1.5i 1.75i 2i
VALUE:\t\c
.br
Test statistic sn is asymptotically Chi-squared with rank(X1) dfs.
The vector of ranks is also returned as component rank.
.DT
See GJKP(1993)
.SH REFERENCES 
       [1] Gutenbrunner,  C.,  J.  Jureckova,   Koenker,  R.  and
       Portnoy,  S.(1993)  "Tests  of Linear Hypotheses  based on
       Regression  Rank   Scores",   Journal   of   Nonparametric
       Statistics, (2), 307-331.

       [2] Koenker, R.W. and d'Orey (1994).  "Remark on  Alg.  AS
       229:  Computing  Dual  Regression Quantiles and Regression
       Rank Scores", Applied Statistics, 43, 410-414.
.SA
rq, rq.ranks
.EX
# Test that covariates 2 and 3 belong in stackloss model using Wilcoxon scores.
rrs.test(stack.x[,1],stack.x[,2:3],stack.loss)
.KW
rank test, quantile regression
.SH Documentation Class
function
.WR
