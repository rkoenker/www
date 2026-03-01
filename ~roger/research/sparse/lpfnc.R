"rqc" <- function(X,y,R,r,tau=.5){
n <- nrow(X)
p <- nrow(X)
m <- nrow(R)
u <- rep(1,n)
a1 <- (1-tau)*u
a2 <- rep(1,m)
b <- t(X)%*%a1
f <- lpfnc(t(X),-y,t(R),-r,b,u,a1,a2)
list(coef = -f$coef, it = f$it)
}
"lpfnc" <-
function(A1,c1,A2,c2,b,u,x1,x2){
#This is a pure R implementation of the inequality constrained interior
#point LP solver:  This was constructed purely to prototype the development
#of fortran versions and should not be taken very seriously for problems of any size.
# NB.  The if(any(is.infinite(q1)))  is tacky and should be fixed in some better way.
	beta <- .9995
	small <- 1e-8
	maxit <- 50
	s <- u-x1
	n1 <- ncol(A1)
	n2 <- ncol(A2)
	y <- coef(lm(c1 ~ t(A1)-1))
	r1 <- c1-t(A1)%*%y
	r2 <- c2-t(A2)%*%y
	z1 <- r1*(r1>0) + small
	w <- z1-r1 + small
	z2 <- rep(1,n2)
	gap <- t(z1)%*%x1 + t(z2)%*%x2 + t(w)%*%s
	it <- 0
	while(gap>small & it<maxit){
		it <- it+1
		q1 <- 1/(z1/x1+w/s)
		#if(any(is.infinite(q1))) 
		#	q1  <- rep(1,length(q1))
		q2 <- x2/z2
		r1 <- c1-t(y%*%A1)
		r2 <- c2-t(y%*%A2)
		r3 <- b-A1%*%x1-A2%*%x2
		Q1 <- diag(as.vector(q1))
		Q2 <- diag(as.vector(q2))
		AQ1 <- A1%*%Q1
		AQ2 <- A2%*%Q2
		AQA <- AQ1%*%t(A1)+AQ2%*%t(A2)
		rhs <- r3+AQ1%*%r1+AQ2%*%r2
		dy <- t(solve(AQA,rhs))
		dx1 <- q1*(t(dy%*%A1)-r1)
		dx2 <- q2*(t(dy%*%A2)-r2)
		ds <- -dx1
		dz1 <- -z1*(1+dx1/x1)
		dz2 <- -z2*(1+dx2/x2)
		dw <- -w*(1+ds/s)
		fx1 <- bound(x1,dx1)
		fx2 <- bound(x2,dx2)
		fz1 <- bound(z1,dz1)
		fz2 <- bound(z2,dz2)
		fs <- bound(s,ds)
		fw <- bound(w,dw)
		fp <- min(pmin(fx1,fs))
		fd <- min(pmin(fw,fz1))
		fp <- min(fp,min(fx2))
		fd <- min(fd,min(fz2))
		fp <- min(min(beta*fp),1)
		fd <- min(min(beta*fd),1)
		if(min(fp,fd)<1){
			mu <- t(z1)%*%x1+t(z2)%*%x2+t(w)%*%s
			g <- t(z1+fd*dz1)%*%(x1+fp*dx1)+ 
				t(z2+fd*dz2)%*%(x2+fp*dx2)+
				t(w+fd*dw)%*%(s+fp*ds)
			mu <- as.vector(mu*(g/mu)^3/(2*n1+n2))
			dxdz1 <- dx1*dz1
			dxdz2 <- dx2*dz2
			dsdw <- ds*dw
			xinv1 <- 1/x1
			xinv2 <- 1/x2
			sinv <- 1/s
			xi1 <- xinv1 * dxdz1 - sinv * dsdw - mu*(xinv1-sinv)
			xi2 <- xinv2 * dxdz2  - mu*xinv2
			rhs <- rhs+ A1%*%(q1*xi1) + A2%*%(q2*xi2)
			dy <- t(solve(AQA,rhs))
			dx1 <- q1*(t(dy%*%A1) - xi1 - r1)
			dx2 <- q2*(t(dy%*%A2) - xi2 - r2)
			ds <- -dx1
			dz1 <- -z1 + xinv1 * (mu - z1*dx1 - dxdz1)
			dz2 <- -z2 + xinv2 * (mu - z2*dx2 - dxdz2)
			dw <- -w + sinv * (mu - w*ds - dsdw)
			fx1 <- bound(x1,dx1)
			fx2 <- bound(x2,dx2)
			fz1 <- bound(z1,dz1)
			fz2 <- bound(z2,dz2)
			fs <- bound(s,ds)
			fw <- bound(w,dw)
			fp <- min(pmin(fx1,fs))
			fd <- min(pmin(fw,fz1))
			fp <- min(fp,min(fx2))
			fd <- min(fd,min(fz2))
			fp <- min(min(beta*fp),1)
			fd <- min(min(beta*fd),1)
			}
		x1 <- x1+fp*dx1
		x2 <- x2+fp*dx2
		z1 <- z1+fd*dz1
		z2 <- z2+fd*dz2
		s <- s+fp*ds
		y <- y+fd*dy
		w <- w+fd*dw
		gap <- t(z1)%*%x1+t(z2)%*%x2+t(w)%*%s
		}
	list(coef = y, it = it)
}
bound <- function(x,dx){
b <- 1e20+0*x
f <- (dx<0)
b[f] <- -x[f]/dx[f]
return(b)
}

