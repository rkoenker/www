# MERDE:  Maximum entropy regularized density estimation 
# Google sez:  "Le mot CRIMERDE n'est pas valide au scrabble."

triweights <- function(x, y = NULL){
    # Riemann weights for scattered data
    n <- nrow(x)
    T <- tri.mesh(x[,1],x[,2])
    T3 <- triangles(T)[, 1:3]
    A <- matrix(0, nrow(T3), n)
    for (i in 1:nrow(T3)) 
        A[i, T3[i, ]] <- abs(det(cbind(x[T3[i, ], ], rep(1, 3))))/2
    apply(A, 2, sum)
}



# See Appendix B of http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
# Or better yet: the original source seems to be Dunavant, D.A.  (1985) High Degree efficient symmetrical
# Gaussian quadrature rules for the triangle, Intl J for Num Methods in Engineering, 21, 1129-1148 

tri.area <- function(v) {
        0.5 * ((v[2, 1] - v[1, 1]) * (v[3, 2] - v[1, 2]) - 
	       (v[3, 1] - v[1, 1]) * (v[2, 2] - v[1, 2]))
    }
GQpts <- function(n){
    switch(n,
	matrix(c(0.33333333333333, 0.33333333333333, 1.00000000000000), 1, 3, byrow = TRUE),
	matrix(c(0.16666666666667, 0.16666666666667, 0.33333333333333,
	    0.16666666666667, 0.66666666666667, 0.33333333333333, 
	    0.66666666666667, 0.16666666666667, 0.33333333333333), 3, 3, byrow = TRUE),
	matrix(c(0.33333333333333, 0.33333333333333, -0.56250000000000,
	    0.20000000000000, 0.20000000000000, 0.52083333333333, 
	    0.20000000000000, 0.60000000000000, 0.52083333333333, 
	    0.60000000000000, 0.20000000000000, 0.52083333333333), 4, 3, byrow = TRUE),
	matrix(c(0.44594849091597, 0.44594849091597, 0.22338158967801,
	    0.44594849091597, 0.10810301816807, 0.22338158967801, 
	    0.10810301816807, 0.44594849091597, 0.22338158967801, 
	    0.09157621350977, 0.09157621350977, 0.10995174365532, 
	    0.09157621350977, 0.81684757298046, 0.10995174365532, 
	    0.81684757298046, 0.09157621350977, 0.10995174365532), 6, 3, byrow = TRUE),
	matrix(c(0.33333333333333, 0.33333333333333, 0.22500000000000,
	    0.47014206410511, 0.47014206410511, 0.13239415278851, 
	    0.47014206410511, 0.05971587178977, 0.13239415278851, 
	    0.05971587178977, 0.47014206410511, 0.13239415278851, 
	    0.10128650732346, 0.10128650732346, 0.12593918054483, 
	    0.10128650732346, 0.79742698535309, 0.12593918054483, 
	    0.79742698535309, 0.10128650732346, 0.12593918054483), 7, 3, byrow = TRUE),
	matrix(c(0.24928674517091, 0.24928674517091, 0.11678627572638,
	    0.24928674517091, 0.50142650965818, 0.11678627572638, 
	    0.50142650965818, 0.24928674517091, 0.11678627572638, 
	    0.06308901449150, 0.06308901449150, 0.05084490637021, 
	    0.06308901449150, 0.87382197101700, 0.05084490637021, 
	    0.87382197101700, 0.06308901449150, 0.05084490637021, 
	    0.31035245103378, 0.63650249912140, 0.08285107561837, 
	    0.63650249912140, 0.05314504984482, 0.08285107561837, 
	    0.05314504984482, 0.31035245103378, 0.08285107561837, 
	    0.63650249912140, 0.31035245103378, 0.08285107561837, 
	    0.31035245103378, 0.05314504984482, 0.08285107561837, 
	    0.05314504984482, 0.63650249912140, 0.08285107561837), 12, 3, byrow = TRUE),
	matrix(c(0.33333333333333, 0.33333333333333, -0.14957004446768,
	    0.26034596607904, 0.26034596607904, 0.17561525743321, 
	    0.26034596607904, 0.47930806784192, 0.17561525743321, 
	    0.47930806784192, 0.26034596607904, 0.17561525743321, 
	    0.06513010290222, 0.06513010290222, 0.05334723560884, 
	    0.06513010290222, 0.86973979419557, 0.05334723560884, 
	    0.86973979419557, 0.06513010290222, 0.05334723560884,
	    0.31286549600487, 0.63844418856981, 0.07711376089026, 
	    0.63844418856981, 0.04869031542532, 0.07711376089026, 
	    0.04869031542532, 0.31286549600487, 0.07711376089026, 
	    0.63844418856981, 0.31286549600487, 0.07711376089026, 
	    0.31286549600487, 0.04869031542532, 0.07711376089026, 
	    0.04869031542532, 0.63844418856981, 0.07711376089026), 13, 3, byrow = TRUE),
	matrix(c(0.33333333333333, 0.33333333333333, 0.14431560767779,
	    0.45929258829272, 0.45929258829272, 0.09509163426728, 
	    0.45929258829272, 0.08141482341455, 0.09509163426728, 
	    0.08141482341455, 0.45929258829272, 0.09509163426728, 
	    0.17056930775176, 0.17056930775176, 0.10321737053472, 
	    0.17056930775176, 0.65886138449648, 0.10321737053472, 
	    0.65886138449648, 0.17056930775176, 0.10321737053472, 
	    0.05054722831703, 0.05054722831703, 0.03245849762320, 
	    0.05054722831703, 0.89890554336594, 0.03245849762320, 
	    0.89890554336594, 0.05054722831703, 0.03245849762320, 
	    0.26311282963464, 0.72849239295540, 0.02723031417443, 
	    0.72849239295540, 0.00839477740996, 0.02723031417443, 
	    0.00839477740996, 0.26311282963464, 0.02723031417443, 
	    0.72849239295540, 0.26311282963464, 0.02723031417443, 
	    0.26311282963464, 0.00839477740996, 0.02723031417443, 
	    0.00839477740996, 0.72849239295540, 0.02723031417443), 16, 3, byrow = TRUE)
	)
}


Bmake <- function(x, N = 5){
    # Gaussian Quadrature weights and constraint 
    X <- GQpts(N)
    M <- nrow(X)
    P <- cbind(1 - X[,1] - X[,2], X[,1], X[,2])
    n <- nrow(x)
    T <- tri.mesh(x[,1],x[,2])
    T3 <- triangles(T)[, 1:3]
    #BB <- matrix(0, M * nrow(T3),n)
    w <- rep(0, M * nrow(T3))
    ia <- ja <- ra <- rep(0, 3 * nrow(T3) * M)
    plot(T)
    for (i in 1:nrow(T3)) {
	xi <- c(x[T3[i,],1])
	yi <- c(x[T3[i,],2])
	xx <- P %*% xi
	yy <- P %*% yi
	points(xx,yy,cex = 0.2, col = 2)
	b <- bary(xx, yy, cbind(xi,yi))
	k <- (i-1)*M
	K <- (i-1)*3*M
	#BB[(k+1):(k+M),T3[i,]] <- b #dense version
	ia[(K+1):((K+3*M))] <- rep((k+1):(k+M),3)
	ja[(K+1):((K+3*M))] <- rep(T3[i,],each = M)
	ra[(K+1):((K+3*M))] <- b$B
	w[(k+1):((k+M))] <- X[,3]*b$Area
    }
    B <- sparseMatrix(i = ia, j = ja, x = ra, dims = c(M*nrow(T3),n))
    list(B = B, w = w)
}

bary <- function(x, y, v){
    # Returns barycentric coordinates  of (x,y) points
    # relative to triangle formed by coordinates in v.
    n <- length(x)
    if(length(y) != n) stop("x,y lengths don't match")
    Area <- tri.area(v)
    B <- matrix(0, n, 3)
    for(i in 1:n){
	V <- rbind(v, c(x[i],y[i]))
        B[i,1] <- tri.area(V[c(4, 2, 3), ])
        B[i,2] <- tri.area(V[c(1, 4, 3), ])
        B[i,3] <- tri.area(V[c(1, 2, 4), ])
    }
    B <- B/Area
    if (any(B < 0 || B > 1)) stop("barycentric snafu")
    list(B = B, Area = Area)
    }



Amake <- function(x, N = 5) {
    n <- nrow(x)
    m <- n-1
    K <- n * m
    ia <- ja <- ra <- rep(0,4*K)
    i <- 0
    for(k in 1:n){ #Seijo-Sen A matrix
	jj <- 0
	for(j in (1:n)[-k]){
	    jj <- jj + 1
	    i <- i + 1
	    ia[i] <- jj + (k-1)*m
	    ja[i] <- k
	    ra[i] <- x[j,1] - x[k,1]
	    ia[K+i] <- jj + (k-1)*m
	    ja[K+i] <- k + n
	    ra[K+i] <- x[j,2] - x[k,2]
	    ia[2*K+i] <- jj + (k-1)*m
	    ja[2*K+i] <- k + 2*n
	    ra[2*K+i] <- 1 
	    ia[3*K+i] <- jj + (k-1)*m
	    ja[3*K+i] <- j  + 2*n
	    ra[3*K+i] <- -1 
	}
    }
    A <- sparseMatrix(i = ia, j = ja, x = ra, dims = c(K,3*n))
    G <- Bmake(x, N = N)
    p <- ncol(A)
    q <- nrow(A)
    m <- nrow(G$B)
    B <- cbind(Matrix(0,m,2*n), G$B, -Diagonal(m))
    A <- cbind(A, Matrix(0,q,m))
    A <- rbind(A,B)
    list(A = A, w = G$w)
}
merde <- function (x, alpha = 1,  m = NULL, w = NULL, rtol = 1e-8, verb = 5, control = NULL) {
    n <- nrow(x)
    p <- 3 * n
    q <- n * (n-1)
    A <- Amake(x)
    C <- A$w
    A <- A$A
    m <- nrow(A) - q
    L <- t(cbind(Matrix(0, n, 2*n), Diagonal(n), Matrix(0,n,m)))
    if(!length(w)) w <- rep(1,n)/n
    wL <- as.vector(L %*% w)
    LC <- c(rep(1e-8, q), rep(0,m)) 
    UC <- c(rep(Inf, q), rep(0,m))
    LX <- rep(-Inf, p + m)
    UX <- rep(Inf, p + m)
    beta <- alpha/(alpha - 1)
    P <- list(sense = "min")
    P$c <- -wL
    P$A <- A
    P$bx <- rbind(LX, UX)
    P$bc <- rbind(LC, UC)
    opro <- matrix(list(), nrow = 5, ncol = m)
    opro[2, ] <- (p+1):(p+m)
    opro[5, ] <- rep(0, m)
    if (alpha == 1) {
        opro[1, ] <- "exp"
        opro[3, ] <- C
        opro[4, ] <- rep(1, m)
	}
    else if(alpha == 0){
        opro[1, ] <- "log"
        opro[3, ] <- C
        opro[4, ] <- rep(1, n)
	}
    else { 
        opro[1, ] <- "pow"
        opro[3, ] <- abs(1-alpha) * C / alpha
        opro[4, ] <- rep(beta, n)
	}
    P$scopt <- list(opro = opro)
    P$dparam$intpnt_nl_tol_rel_gap <- rtol
    if (length(control)) {
        P$iparam <- control$iparam
        P$dparam <- control$dparam
        P$sparam <- control$sparam
    }
    z <- Rmosek::mosek(P, opts = list(verbose = verb))
    status = z$sol$itr$solsta
    if (status != "OPTIMAL") 
        warning(paste("Solution status = ", status))
    y <- z$sol$itr$xx
    print(range((A %*% y)[1:q]))
    y <- y[(2*n+1):p]
    z <- list(x = x, y = y, status = status)
    class(z) <- "merde"
    z
}

require(Rmosek)
require(tripack)
require(deldir)
require(rgl)
set.seed(14)
x <- matrix(rnorm(100),50,2)
f <- merde(x)
dxyz <- deldir(f$x[,1],f$x[,2],z = f$y)
persp3d(dxyz,smooth = FALSE)

