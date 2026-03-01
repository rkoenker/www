"qss"<-
function(x, y, w = rep(1, length(y)), theta = .5, lambda = -1, method = "bc", penalty="Loo", constraint="none",
	coef = rep(0, nvar), pswitch = F, maxiter = 500, penpar = 1/lam^.5, 
	tmin = 1/(ny + nh), kmax = sqrt(ny))
{
#Compute linear quantile smoothing spline with L1 penalty and
#penalty parameter lambda
#Solution is a linear spline: 
#	g(x) = a_i * (x - x_i) + b_i 
#Inputs:
#	x--explanatory variable(s)
#	y--explained variable
#	w--weights
#	theta--desired quantile
#	lambda--penalty parameter
#Outputs:
#	coef--smoothed value at each design point with coef[i] <-  b_i
#       resid--vector of residuals
#	exit--exit code for "br"
#       ifl--exit code for "bc"
#       niter--number of iteration taken for "bc"
#       fidel--fidelity
#       rough--roughness
#       k--dimension
#       obj--value of the objective function
#	sol--quantile smoothing spline process
#
	big <- .Machine$double.xmax
	single.big <- .Machine$single.xmax
	eps <- .Machine$double.eps
	zero <- 0
	alll <- F
	allt <- F
	if(constraint != "none" & constraint != "increase" & constraint != "decrease" & constraint != "convex" & constraint != "concave")
		stop("constraint has to be none, increase, decrease, convex or concave")
	if(penalty != "L1" & penalty != "Loo")
		stop("penalty has to be either L1 or Loo")
	if(method != "br" & method != "bc")
		stop("method has to be either br or bc")
	if(penalty=="Loo" & method=="br")
		stop("Use bc method for Loo penalty")
	if(constraint != "none" & method != "bc")
		stop("Use bc method for constrained problems")
	if(is.matrix(x)){
		if(constraint!="none")
			stop("Can't impose additional constraints on additive splines yet")
		k <- ncol(x)
		if(!(nrow(x) == length(y) & nrow(x) == length(w)))
			stop("x, y and w do not have compatible dimensions")
		if(length(lambda) != k)
			stop("lambda doesn't have the same length as the column of x")
		if(min(lambda)<0)
			stop("no lambda can be negative")
		if(theta < 0 | theta > 1)
			stop("Can't perform parametric programming for additive splines yet")
		}
	else{
		if(!(length(x) == length(y) & length(x) == length(w)))
			stop("x, y and w do not have equal length")
		if(lambda < 0) {
			alll <- T
			penpar <- 1/eps
			}
		else alll <- F
		if(theta < 0 | theta > 1){
			allt <- T
			}
		else allt <- F
		if(allt & alll)
			stop("Cant' solve for all t and lambda simultaneously")
		}
	ny <- length(y)
	if(!is.matrix(x) & constraint=="none" & penalty=="L1"){
		ox <- order(x)
		xun <- unique(x[ox])
		h <- diff(xun)
		nh <- length(h)
		B <- diag(1/h,nrow=nh)
		B <- cbind(0,B)
		B[row(B) == col(B)] <- -1/h
		D <- diag(1,nrow=(nh-1))
		D <- cbind(0,D)
		D[row(D)==col(D)] <- -1
		A <- D %*% B
		X <- matrix(0, ny, nh + 1)
		X[cbind(1:ny, category(x))] <- 1	#fidelity part of design matrix
		X <- rbind(diag(w) %*% X, A)	#whole X matrix w/o lambda
		Y <- c(w * y, rep(0, nh-1))
		nobs <- nrow(X)
		nvar <- ncol(X)
		nrq <- ny
		nl1 <- nh-1
		neqc <- zero
		niqc <- zero
		lam <- lambda
		ncoef <- nvar
		}
	if(!is.matrix(x) & constraint=="none" & penalty=="Loo"){
		ox <- order(x)
		xun <- unique(x[ox])
		h <- diff(xun)
		nh <- length(h)
		D <- diag(h, nrow = nh)
	        D[row(D) == col(D) + 1] <- h[1:(nh - 1)]
       		D[1, 1] <- 1
	       	B <- diag(1/h, nrow = nh)
       		B[row(B) == col(B) + 1] <-  - (1/h[1:(nh - 1)] + 1/h[2:(nh)])
	       	B[row(B) == col(B) + 2] <- 1/h[2:(nh - 1)]
       		B <- cbind(c(0, 1/h[1], rep(0, nh - 2)), B)
	       	B[1,  ] <- 0
	       	B <- cbind(0, B)
	       	B[1, 1] <- 1
	       	A <- solve(D) %*% B
       	 	X <- matrix(0, ny, nh + 1)
	       	X[cbind(1:ny, category(x))] <- 1        #fidelity part of design matrix
	       	X <- rbind(cbind(0, diag(w) %*% X, 0),c(rep(0,(nh+2)),1))
	       	Y <- c(w * y,0,rep(0,2*nh))
	       	X <- rbind(X, cbind(A, 1), cbind(- A, 1))
		nobs <- nrow(X)
		nvar <- ncol(X)
		nrq <- ny
		nl1 <- 1
		neqc <- zero
		niqc <- 2*nrow(A)
		lam <- lambda
		ncoef <- nvar-1
		}
	if(!is.matrix(x) & constraint=="increase" & penalty=="L1"){
		ox <- order(x)
	        xun <- unique(x[ox])
	        h <- diff(xun)
	        nh <- len(h)
		B <- diag(1/h,nrow=nh)
                B <- cbind(0,B)
                B[row(B) == col(B)] <- -1/h
                D <- diag(1,nrow=(nh-1))
                D <- cbind(0,D)
                D[row(D)==col(D)] <- -1
                A <- D %*% B
		G <- diag(1/h)
		G <- cbind(0,G)
		G[row(G) == col(G)] <- -1/h
		X <- matrix(0, ny, nh + 1)
                X[cbind(1:ny, category(x))] <- 1        #fidelity part of design matrix
                X <- rbind(diag(w) %*% X, A, G)    #whole X matrix w/o lambda
                Y <- c(w * y, rep(0, nh-1), rep(0, nh))
		nobs <- nrow(X)
		nvar <- ncol(X)
		nrq <- ny
		nl1 <- nh-1
		neqc <- zero
		niqc <- nh 
		lam <- lambda
		ncoef <- nvar
		}
	if(!is.matrix(x) & constraint=="decrease" & penalty=="L1"){
		ox <- order(x)
	        xun <- unique(x[ox])
	        h <- diff(xun)
	        nh <- len(h)
		B <- diag(1/h,nrow=nh)
                B <- cbind(0,B)
                B[row(B) == col(B)] <- -1/h
                D <- diag(1,nrow=(nh-1))
                D <- cbind(0,D)
                D[row(D)==col(D)] <- -1
                A <- D %*% B
		G <- diag(1/h)
		G <- cbind(0,G)
		G[row(G) == col(G)] <- -1/h
		X <- matrix(0, ny, nh + 1)
                X[cbind(1:ny, category(x))] <- 1        #fidelity part of design matrix
                X <- rbind(diag(w) %*% X, A, -G)    #whole X matrix w/o lambda
                Y <- c(w * y, rep(0, nh-1), rep(0, nh))
		nobs <- nrow(X)
		nvar <- ncol(X)
		nrq <- ny
		nl1 <- nh-1
		neqc <- zero
		niqc <- nh 
		lam <- lambda
		ncoef <- nvar
		}
	if(!is.matrix(x) & constraint=="convex" & penalty=="L1"){
		ox <- order(x)
		xun <- unique(x[ox])
		h <- diff(xun)
		nh <- length(h)
		B <- diag(1/h,nrow=nh)
		B <- cbind(0,B)
		B[row(B) == col(B)] <- -1/h
		D <- diag(1,nrow=(nh-1))
		D <- cbind(0,D)
		D[row(D)==col(D)] <- -1
		A <- D %*% B
		X <- matrix(0, ny, nh + 1)
		X[cbind(1:ny, category(x))] <- 1	#fidelity part of design matrix
		X <- rbind(diag(w) %*% X, A, A)	#whole X matrix w/o lambda
		Y <- c(w * y, rep(0, nh-1))
		nobs <- nrow(X)
		nvar <- ncol(X)
		nrq <- ny
		nl1 <- nh-1
		neqc <- zero
		niqc <- nh-1
		lam <- lambda
		ncoef <- nvar
		}
	if(!is.matrix(x) & constraint=="concave" & penalty=="L1"){
		ox <- order(x)
		xun <- unique(x[ox])
		h <- diff(xun)
		nh <- length(h)
		B <- diag(1/h,nrow=nh)
		B <- cbind(0,B)
		B[row(B) == col(B)] <- -1/h
		D <- diag(1,nrow=(nh-1))
		D <- cbind(0,D)
		D[row(D)==col(D)] <- -1
		A <- D %*% B
		X <- matrix(0, ny, nh + 1)
		X[cbind(1:ny, category(x))] <- 1	#fidelity part of design matrix
		X <- rbind(diag(w) %*% X, A, -A)	#whole X matrix w/o lambda
		Y <- c(w * y, rep(0, nh-1))
		nobs <- nrow(X)
		nvar <- ncol(X)
		nrq <- ny
		nl1 <- nh-1
		neqc <- zero
		niqc <- nh-1
		lam <- lambda
		ncoef <- nvar
		}
	if(!is.matrix(x) & constraint=="increase" & penalty=="Loo"){
		ox <- order(x)
		xun <- unique(x[ox])
		h <- diff(xun)
		nh <- length(h)
		D <- diag(h, nrow = nh)
	        D[row(D) == col(D) + 1] <- h[1:(nh - 1)]
       		D[1, 1] <- 1
	        B <- diag(1/h, nrow = nh)
	        B[row(B) == col(B) + 1] <-  - (1/h[1:(nh - 1)] + 1/h[2:(nh)])
	        B[row(B) == col(B) + 2] <- 1/h[2:(nh - 1)]
	        B <- cbind(c(0, 1/h[1], rep(0, nh - 2)), B)
	        B[1,  ] <- 0
	        B <- cbind(0, B)
	        B[1, 1] <- 1
	        A <- solve(D) %*% B
	        G <- diag(-1/h)
	        G <- cbind(G, 0)
	        G[row(G) + 1 == col(G)] <- 1/h
	        G <- cbind(0, G)
	        G <- rbind(G, G[1,  ])
	        C <- diag(h, nrow = nh)
	        C <- rbind(C, c( - h[1], -2 * h[2:nh]))
	        G <- G - C %*% A
	        X <- matrix(0, ny, nh + 1)
	        X[cbind(1:ny, category(x))] <- 1     #fidelity part of design matrix
	        X <- rbind(cbind(0, diag(w) %*% X, 0),c(rep(0,(nh+2)),1))
	        Y <- c(w * y,0)
	        X <- rbind(X, cbind(A, 1), cbind(- A, 1), cbind(G, 0))
		Y <- c(Y,rep(0,2*nh),rep(0,nh+1))
		nobs <- nrow(X)
		nvar <- ncol(X)
		nrq <- ny
		nl1 <- 1
		neqc <- zero
		niqc <- 2 * nrow(A) + nrow(G)
		lam <- lambda
		ncoef <- nvar - 1
		}
	if(!is.matrix(x) & constraint=="decrease" & penalty=="Loo"){
		ox <- order(x)
		xun <- unique(x[ox])
		h <- diff(xun)
		nh <- length(h)
		D <- diag(h, nrow = nh)
	        D[row(D) == col(D) + 1] <- h[1:(nh - 1)]
       		D[1, 1] <- 1
	        B <- diag(1/h, nrow = nh)
	        B[row(B) == col(B) + 1] <-  - (1/h[1:(nh - 1)] + 1/h[2:(nh)])
	        B[row(B) == col(B) + 2] <- 1/h[2:(nh - 1)]
	        B <- cbind(c(0, 1/h[1], rep(0, nh - 2)), B)
	        B[1,  ] <- 0
	        B <- cbind(0, B)
	        B[1, 1] <- 1
	        A <- solve(D) %*% B
	        G <- diag(-1/h)
	        G <- cbind(G, 0)
	        G[row(G) + 1 == col(G)] <- 1/h
	        G <- cbind(0, G)
	        G <- rbind(G, G[1,  ])
	        C <- diag(h, nrow = nh)
	        C <- rbind(C, c( - h[1], -2 * h[2:nh]))
	        G <- G - C %*% A
	        X <- matrix(0, ny, nh + 1)
	        X[cbind(1:ny, category(x))] <- 1     #fidelity part of design matrix
	        X <- rbind(cbind(0, diag(w) %*% X, 0),c(rep(0,(nh+2)),1))
	        Y <- c(w * y,0)
	        X <- rbind(X, cbind(A, 1), cbind(- A, 1), -cbind(G, 0))
		Y <- c(Y,rep(0,2*nh),rep(0,nh+1))
		nobs <- nrow(X)
		nvar <- ncol(X)
		nrq <- ny
		nl1 <- 1
		neqc <- zero
		niqc <- 2 * nrow(A) + nrow(G)
		lam <- lambda
		ncoef <- nvar - 1
		}
	if(!is.matrix(x) & constraint=="convex" & penalty=="Loo"){
		ox <- order(x)
		xun <- unique(x[ox])
		h <- diff(xun)
		nh <- length(h)
		D <- diag(h, nrow = nh)
	        D[row(D) == col(D) + 1] <- h[1:(nh - 1)]
       		D[1, 1] <- 1
	       	B <- diag(1/h, nrow = nh)
       		B[row(B) == col(B) + 1] <-  - (1/h[1:(nh - 1)] + 1/h[2:(nh)])
	       	B[row(B) == col(B) + 2] <- 1/h[2:(nh - 1)]
       		B <- cbind(c(0, 1/h[1], rep(0, nh - 2)), B)
	       	B[1,  ] <- 0
	       	B <- cbind(0, B)
	       	B[1, 1] <- 1
	       	A <- solve(D) %*% B
       	 	X <- matrix(0, ny, nh + 1)
	       	X[cbind(1:ny, category(x))] <- 1        #fidelity part of design matrix
	       	X <- rbind(cbind(0, diag(w) %*% X, 0),c(rep(0,(nh+2)),1))
	       	Y <- c(w * y,0,rep(0,2*nh))
	       	X <- rbind(X, cbind(A, 1), cbind(- A, 1), cbind(A, 0))
		nobs <- nrow(X)
		nvar <- ncol(X)
		nrq <- ny
		nl1 <- 1
		neqc <- zero
		niqc <- 3*nrow(A)
		lam <- lambda
		ncoef <- nvar-1
		}
	if(!is.matrix(x) & constraint=="concave" & penalty=="Loo"){
		ox <- order(x)
		xun <- unique(x[ox])
		h <- diff(xun)
		nh <- length(h)
		D <- diag(h, nrow = nh)
	        D[row(D) == col(D) + 1] <- h[1:(nh - 1)]
       		D[1, 1] <- 1
	       	B <- diag(1/h, nrow = nh)
       		B[row(B) == col(B) + 1] <-  - (1/h[1:(nh - 1)] + 1/h[2:(nh)])
	       	B[row(B) == col(B) + 2] <- 1/h[2:(nh - 1)]
       		B <- cbind(c(0, 1/h[1], rep(0, nh - 2)), B)
	       	B[1,  ] <- 0
	       	B <- cbind(0, B)
	       	B[1, 1] <- 1
	       	A <- solve(D) %*% B
       	 	X <- matrix(0, ny, nh + 1)
	       	X[cbind(1:ny, category(x))] <- 1        #fidelity part of design matrix
	       	X <- rbind(cbind(0, diag(w) %*% X, 0),c(rep(0,(nh+2)),1))
	       	Y <- c(w * y,0,rep(0,2*nh))
	       	X <- rbind(X, cbind(A, 1), cbind(- A, 1), -cbind(A, 0))
		nobs <- nrow(X)
		nvar <- ncol(X)
		nrq <- ny
		nl1 <- 1
		neqc <- zero
		niqc <- 3*nrow(A)
		lam <- lambda
		ncoef <- nvar-1
		}
	if(is.matrix(x) & constraint=="none" & penalty=="L1"){
		for(i in 1:k){
			ox <- order(x[,i])
			xun <- unique((x[,i])[ox])
			h <- diff(xun)
			nh <- len(h)
			B <- diag(1/h,nrow=nh)
	       		B <- cbind(0,B)
	       		B[row(B) == col(B)] <- -1/h
       			D <- diag(1,nrow=(nh-1))
       			D <- cbind(0,D)
		       	D[row(D)==col(D)] <- -1
		       	A <- D %*% B
		       	X <- matrix(0, ny, nh + 1)
		       	X[cbind(1:ny, category(x[,i]))] <- 1    #fidelity part of design matrix
			if(i==1){
				Xtilde <- rbind(diag(w) %*% X,lambda[i]*A)
				Xconst <- matrix(c(1,rep(0,nh)),nrow=1)
				next
				}
			Xtilde <- cbind(rbind(Xtilde,matrix(0,nrow=nh-1,ncol=ncol(Xtilde))),
				rbind(diag(w) %*% X,matrix(0,nrow=(nrow(Xtilde)-ny),
				ncol=nh+1),lambda[i]*A)) #whole X  matrix with lambda
			Xconst <- rbind(cbind(Xconst,matrix(0,nrow=i-1,ncol=nh+1)),c(rep(0,ncol(Xconst)),c(1,rep(0,nh))))
			}
		nobs <- nrow(Xtilde)
		X <- rbind(cbind(Xtilde,c(rep(1,ny),rep(0,nobs-ny))),cbind(Xconst,0))
		nobs <- nrow(X)
		nvar <- ncol(X)
		Y <- c(w * y, rep(0, nobs-ny))
		nrq <- ny
		nl1 <- nobs-ny-k
		neqc <- k
		niqc <- zero
		penpar <- 1/median(lambda)
		lam <- 1
		ncoef <- nvar
		}
	if(is.matrix(x) & constraint=="none" & penalty=="Loo"){
		for(i in 1:k){
			ox <- order(x[,i])
			xun <- unique((x[,i])[ox])
			h <- diff(xun)
			nh <- len(h)
			D <- diag(h, nrow = nh)
			D[row(D) == col(D) + 1] <- h[1:(nh - 1)]
			D[1, 1] <- 1
			B <- diag(1/h, nrow = nh)
			B[row(B) == col(B) + 1] <-  - (1/h[1:(nh - 1)] + 1/h[2:(nh)])
			B[row(B) == col(B) + 2] <- 1/h[2:(nh - 1)]
			B <- cbind(c(0, 1/h[1], rep(0, nh - 2)), B)
			B[1,  ] <- 0
			B <- cbind(0, B)
			B[1, 1] <- 1
			A <- solve(D) %*% B
			X <- matrix(0, ny, nh + 1)
			X[cbind(1:ny, category(x[,i]))] <- 1	#fidelity part of design matrix
			if(i==1){
				Xtilde <- rbind(cbind(0, diag(w) %*% X, 0),
				c(rep(0,nh+2), lambda[i]), cbind(A,1), cbind(-A,1))
				Xconst <- matrix(c(0,1,rep(0,nh+1)),nrow=1)
				next
				}
			Xtilde <- cbind(rbind(Xtilde,matrix(0,nrow=2*nh,ncol=ncol(Xtilde))),
			rbind(cbind(0,diag(w) %*% X, 0), c(rep(0, nh+2), 
			lambda[i]), matrix(0,nrow=(nrow(Xtilde)-(ny+1)),
			ncol=nh+3),cbind(A,1), cbind(-A,1))) #whole X  matrix with lambda
			Xconst <- rbind(cbind(Xconst,matrix(0,nrow=i-1,ncol=nh+3)),c(rep(0,ncol(Xconst)),c(0,1,rep(0,nh+1))))
			}
		nobs <- nrow(Xtilde)
		X <- cbind(rbind(Xtilde[1:(ny+1),],Xconst,Xtilde[(ny+2):nobs,]),c(rep(1,ny),rep(0,nobs+k-ny)))
		nobs <- nrow(X)
		nvar <- ncol(X)
		Y <- c(w * y, rep(0, nobs-ny))
		nrq <- ny
		nl1 <- 1
		neqc <- k
		niqc <- nobs-ny-1-k
		penpar <- 1/median(lambda)
		lam <- 1
		ncoef <- nvar
		}
	storage.mode(X) <- "single"
	d <- matrix(0, nobs + 5, nvar + 2)
	storage.mode(d) <- "double"
	sd <- matrix(0, nobs + 5, nvar + 2)
	storage.mode(sd) <- "single"
	nj0 <- 1
	if(allt|alll)
		nj0 <- 10*nobs
	sol <- matrix(0, nvar + 6, nj0)
	dimnames(sol)<-list(c("quantile","lambda","exit code","fidelity","penalty","dimension",rep("coef",nvar)),NULL)
	toler <- 1e-08	#Note: toler might need fine tuning; can't be too
#small for tnt > t; too big will miss some t
	if(method == "br") {
		storage.mode(sol) <- "single"
		z <- .Fortran("drqssbr",
			as.integer(nobs),
			as.integer(nvar),
			as.integer(nobs + 5),
			as.integer(nvar + 2),
			X = X,
			as.single(Y),
			as.single(theta),
			as.single(toler),
			coef = single(nvar),
			resid = single(nobs),
			integer(nobs),
			d,
			double(nobs),
			as.integer(nj0),
			sol = sol,
			nj1 = integer(1),
			as.single(lam),
			sd = sd,
			as.single(single.big),
			as.single(eps),
			k = integer(1),	
			as.integer(kmax))
		if(allt | alll) {
			z$sol <- z$sol[, 1:z$nj1]
			return(sol = z$sol)
		}
		else {
			return(coef = z$coef[1:ncoef], 
				resid = z$resid[1:nrq], 
				exit = z$sd[nobs + 2,nvar + 1], 
				fidel = sum((theta-(z$resid[1:nrq]<0))*(z$resid[1:nrq]))*2,
				rough = sum(abs(z$resid[(nrq+1):(nrq+nl1)])),
				obj = z$sd[nobs + 1, nvar + 1], 
				k = z$k,
				penalty = penalty)
		}
	}
	if(method=="bc"){
		storage.mode(sol) <- "double"
		z <- .Fortran("drqssbc",
			as.integer(c(nrq, nl1, neqc, niqc, nvar)),
			integer(1),
			ifl = integer(1),
			as.integer(maxiter),
			as.logical(pswitch),
			X = as.double(t(X)),
			as.integer(nvar),
			coef = as.double(coef),
			as.double(Y),
			obj = double(1),
			resid = double(nobs),
			integer(nobs),
			double((3 * nvar * nvar + 13 * nvar + 2)/2 + 2 * nobs),

				nj1 = integer(1),
			as.integer(nj0),
			sol = sol,
			as.double(c(theta, lam)),
			as.double(toler),
			as.double(big),
			as.double(eps),
			as.double(penpar),
			icyc = integer(1),
			as.double(tmin),
			k = integer(1),
			as.integer(kmax))
		if(allt | alll) {
			storage.mode(z$sol) <- "single"
			return(sol = z$sol[, 1:z$nj1])
		}
		else {
			return(
				coef = as.single(z$coef[1:ncoef]), 
				resid = as.single(-z$resid[1:nrq]), 
				ifl = z$ifl, 
				niter = z$icyc, 
				fidel= as.single(sum((theta-(-z$resid[1:nrq]<0))*(-z$resid[1:nrq]))*2),
				rough= as.single(sum(abs(z$resid[(nrq+1):(nrq+nl1)])))/lam,
				obj = as.single(z$obj), 
				k = z$k,
				penalty = penalty)
			}
	}
}
