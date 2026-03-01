.BG
.FN nlrq
.TL
function to compute nonlinear quantile regression estimates
.DN
This function implements an S+ version of an interior point method
for computing the solution to quantile regression problems which
are nonlinear in the parameters. 
.CS
nlrq(model, t, k=2, theta=0.5, nit.max=100, eps=1e-07, beta=0.97)
.RA
.AG model
is a S function which maps parameters to (residual,Jacobian) space
i.e. given a t, it returns a pair m$f and m$J which are the vector
of residuals and the Jacobian matrix of the model at this t.
.AG t
the initial value of the parameter
.OA
.AG k
the number of iterations of the Meketon algorithm to be done in each
step, usually 2 is reasonable, occasionally it may be helpful to set k=1
.AG theta
The quantile being estimated
.AG nit.max
maximum number of allowed iterations
.AG eps
tolerance for convergence of the algorithm
.AG beta
a shrinkage parameter which controls the recentering process in the
interior point algorithm.
.RT
The function returns three components:  coef-the optimal value of the
parameter vector, obj-the optimal value of the objective function, and
nit- the number of iterations taken
.SH REFERENCES
Koenker, R. and Park, BJ(1994) An Interior Point Algorithm for
Nonlinear Quantile Regression, Journal of Econometrics, forthcoming.
.SA
rq
.EX
brown<- function(b)
{
	t <- (1:20)/5
	f <- (b[1] + t * b[2] - exp(t))^2 + (b[3] + 
		b[4] * sin(t) - cos(t))^2
	J <- matrix(0, 20, 4)
	J[, 1] <- -2 * (b[1] + t * b[2] - exp(t))
	J[, 2] <- -2 * t * (b[1] + t * b[2] - exp(t))
	J[, 3] <- -2 * (b[3] + b[4] * sin(t) - cos(t))
	J[, 4] <- -2 * sin(t) * (b[3] + b[4] * sin(t) - cos(t))
	return(f, J)
}
nlrq(brown,c(25,5,-5,-1))
.WR
