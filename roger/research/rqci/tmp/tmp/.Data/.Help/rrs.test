.BG
.FN rrs.test
.TL
Function to compute regression rankscore test of linear hypothesis
.DN
Tests the hypothesis that b_1 = 0 in the quantile regression model

	y = X_0 b_0 + X_1 b_1 + u 
.CS
rrs.test(x0, x1, y, v, score="wilcoxon")
.RA
.AG x0
matrix of maintained regressors, a column of ones is appended automatically
.AG x1
matrix of regressors under test
.OA
.AG y
response variable may be omitted if v is provided
.AG v
regression quantile structure from rq(x0,y)
.AG score
score function for test (see ranks())
.RT
.RC sn
Test statistic is asymptotically Chi-squared with rank(X1) degrees of freedom.
.RC rank
vector of ranks
.SE
.DT
.SH REFERENCES
[1] Gutenbrunner, C., J. Jureckova,  Koenker, R. and  Portnoy, S.(1993)
"Tests of Linear Hypotheses  based on Regression Rank Scores",
Journal of Nonparametric Statistics, (2), 307-331.

[2] Koenker, R.W. and d'Orey (1994).  "Remark on Alg. AS 229: Computing Dual
Regression Quantiles and Regression Rank Scores", Applied Statistics, 43, 410-414.

.SA
rq, ranks
.EX
rrs.test(stack.x[,1:2],stack.x[,3],stack.loss)
.KW regression quantiles, ranks, regresion rankscores
.WR
