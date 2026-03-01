# Example Plot for simulation a la Ruppert Carroll and Wand Section 17.5.1

# The nonparametric target function(s) as j gets larger g0 has larger total variation
g0 <- function(x, j) sqrt(x *(1-x)) * sin(2*pi*(1 + 2^((9-4*j)/5))/(x + 2^((9-4*j)/5)))



fitrqss <- function(x,y,g0 = NULL, title = NULL){
     g <- function(lam,y,x) AIC(rqss(y ~ qss(x, lambda = lam)),k = -1) 
     lamstar <- optimize(g, interval = c(0.001, .5), x = x, y = y)
     f <- rqss(y ~ qss(x, lambda = lamstar$min))
     plot(x,y,xlab = "x", ylab = "", type = "n")
     B <- plot(f, add = TRUE, bands = "both",title = title)
     bandg <- B[[1]]
     gtrue <- g0(bandg$x,4)
     lines(bandg$x, gtrue, col = "red")
     points(x,y,col = "blue", cex = .25)
     ghat <- (bandg$bhi[,1] + bandg$blo[,1])/2
     lines(bandg$x, ghat, col = "black")
     }
require(quantreg)

n <- 400
x <- sort(runif(n))
sigma <- 0.2
options(warn = -1) # turn off warnings about tiny diagonals
y <- g0(x,4) + rnorm(n)*sigma
fitrqss(x,y, g0, title = "Median Estimate")
