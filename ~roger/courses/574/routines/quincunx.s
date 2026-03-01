function(balls = 100, layers = 30, speed = 4, cex = 1)
{
	bins <- NULL
	plot(c( - layers, layers), c(-2, layers), type = "n", bty = "n", xaxt
		 = "n", yaxt = "n", ylab = "", xlab = "")
	title("Quincunx")
	x <-  - layers:layers
	y <- abs(x)
	b <- 2 * (1:(layers/2)) - 1
	b <- c( - b, b)
	segments(b, rep(-2, length(b)), b, rep(0, length(b)))
	for(i in 1:layers) {
		points(( - (i - 1):(i - 1))[(1:i) * 2], rep(layers - i, i))
	}
	for(i in 1:(balls/speed)) {
		path <- matrix(2 * (runif(speed * layers) > 0.5) - 1, layers, 
			speed)
		x <- apply(path, 2, "cumsum")
		old.y <- rep(layers - 1.5, speed)
		old.x <- rep(0, speed)
		for(k in 1:layers) {
			eps <- rep(rnorm(1)/10, speed)
			points(x[k,  ] + eps, rep(layers - k - 1.5, speed), cex
				 = cex, pch = "o")
			old.y <- rep(layers - k - 1.5, speed)
			old.x <- x[k,  ] + eps
			points(old.x, old.y, cex = cex, pch = "o", col = 0)
		}
		hgt <- i/(balls/speed)
		points(x[layers,  ] + runif(speed)/1.5, -1.5 + (hgt * runif(
			speed)), cex = cex, pch = "o")
		bins <- c(bins, x[layers,  ])
	}
	return(table(bins))
}
