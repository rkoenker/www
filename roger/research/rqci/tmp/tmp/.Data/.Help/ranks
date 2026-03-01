.BG
.FN ranks
.TL
Function to compute ranks from the dual (regression rankscore) process
.DN
.CS
ranks(v, score="wilcoxon")
.RA
.AG v
regression quantile structure for the model of interest
.OA
.AG score
The score function desired.  Currently implemented score functions are
Wilcoxon, Normal, and Sign which are asymptotically optimal for the
logistic, Gaussian and Laplace error models respectively.
.RT
The function returns two components one is the ranks, the other is a
scale factor which is the L_2 norm of the score function.  All score
functions should be normalized to have mean zero.
.SE
.DT
.SH REFERENCES
Gutenbrunner, C., J. Jureckova,  Koenker, R. and  Portnoy, S.(1993)
Tests of Linear Hypotheses  based on Regression Rank Scores",
Journal of Nonparametric Statistics, (2), 307-331.
.SA
See also rq, rrs.test
.EX
ranks(rq(stack.x,stack.loss))
.KW regression quantiles, ranks, regression rankscores
.WR
