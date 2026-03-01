"qss.fit"<-
function(x, g, z)
{
#compute fitted value of the quantile smoothing splines at points z 
#given knots at x
#Inputs:
#	qss.o--output from qss
#	z--vector or list of points where smoothed values are desired 
#Output:
#	returns fitted values at z
#
	penalty <- qss.o$penalty
	BIG <- .Machine$single.xmax
	if(penalty=="Loo"){
		if(!is.matrix(x)){
			xo_order(x)
			x_unique(x[xo])
			n_length(x)
			h <- c(0, diff(x))
			a <- rep(0, n + 1)
			a[2] <- g[1]
			g[1] <- g[2]
			b <- rep( - a[2] * h[2] + (g[3] - g[2])/h[2], n + 1)
			for(i in 3:(n + 1)) {
				b[i] <- (2 * (g[i] - g[i - 1]))/h[i - 1] - b[i - 1]
				}
			a[3:n] <- diff(b[3:(n + 1)])/(2 * h[3:n])
			k <- cut(z, c(x, BIG)) + 1	#obtain bin numbers for z
			k[z <= x[1]] <- 1
			x <- c(x[1], x)
			dz <- z - x[k]
			return(g[k] + b[k] * dz + a[k] * dz^2)
			}
		else if(!is.list(z))
			stop("z should have as many components as ncol(x)")
		else if(ncol(x) > 2)
			stop("Can't perform slicing for x having more than two independent variables yet")
		else {
			x1.uni <- unique((x[,1])[order(x[,1])])
			x2.uni <- unique((x[,2])[order(x[,2])])
			ncoef <- length(g)
			delta <- g[ncoef]
			zz <- matrix(0,nrow=length(z[[1]]),ncol=length(z[[2]]))
			z1.o <- qss.fit(x1.uni,g[1:(length(x1.uni)+1)],z[[1]])
			z2.o <- qss.fit(x2.uni,g[(length(x1.uni)+3):(length(x1.uni)+3+length(x2.uni))],z[[2]])
			for(i in 1:length(z[[1]])){
				for(j in 1:length(z[[2]])){
					zz[i,j] <- delta + z1.o[i] + z2.o[j]
					}
				}
			return(zz)
			}
		}
	if(penalty=="L1"){
		if(!is.matrix(x)){
			xo_order(x)
			x_unique(x[xo])
			n_length(x)
			h <-  diff(x)
			nh <- length(h)
			b <- diff(g)/h
			g <- c(g[1],g)
			b <- c(b[1],b,b[nh])
			k <- cut(z, c(x, BIG)) + 1	#obtain bin numbers for z
			k[z <= x[1]] <- 1
			x <- c(x[1], x)
			dz <- z - x[k]
			return(g[k] + b[k] * dz)
			}
		else if(!is.list(z))
			stop("z should be a list with $z1 and $z2")
		else {
			x1.uni <- unique((x[,1])[order(x[,1])])
			x2.uni <- unique((x[,2])[order(x[,2])])
			ncoef <- length(g)
			delta <- g[ncoef]
			zz <- matrix(0,nrow=length(z[[1]]),ncol=length(z[[2]]))
			z1.o <- qss.fit(x1.uni,g[1:length(x1.uni)],z[[1]])
			z2.o <- qss.fit(x2.uni,g[(length(x1.uni)+1):(length(x1.uni)+length(x2.uni))],z[[2]])
			for(i in 1:length(z[[1]])){
				for(j in 1:length(z[[2]])){
					zz[i,j] <- delta + z1.o[i] + z2.o[j]
					}
				}
			return(zz)
			}
		}
}
