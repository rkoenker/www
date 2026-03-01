.BG
.FN trq.print
.TL
Print a Trimmed Regression Quantile Summary
.DN
Prints out a summary of a trimmed quantile regression (using trq) including the standard error of the residuals, the overall F-statistic and t-statistics and standard errors for each coefficient.
.CS
trq.print(trq.out, digits=4)
.RA
.AG trq.out
a list like the output from trq.
.OA
.AG digits
the number of significant digits that should be printed.
.RT
an invisible list with the following components:
.RC summary
a vector containig the winsorized standard error of regression, 
the number of observations, the overall F-statistic for the regression, 
the degrees of freedom and the p-value for the F-statistic. 
When a1+a2=1, residual winsorized standard
error is the the square root of a Siddiqui-type estimate of the asymtotic
variance of the a1-quantile.
.RC coef.table
a matrix with columns containing the regrssion coefficients, their standard error,the t-statistic for testing if the coefficient os zero and the two sided p-value for the t statistic.
.SE
a summary of the regression (like the returned value) is printed out.
.SA
trq, lsfit, and ls.print for further reference.
.EX
trq.print(trq(x,y))
.KW ~keyword
.WR
