# Very Cheezy attempt to reproduce Fig 1.1 of Parikh and Boyd "Proximal Algorithms"
Pf <- function(v,lam = 3.8){
    g <- function(x,v, lambda) sum(x^2) + lambda * sum((x - v)^2)
    optim(c(2,2), g, v = v, lambda = lam)$par
}
x <- 1:100/20
y <- 1:100/20
pdf("Egprox.pdf", height = 7, width = 7)
plot(x, y, type = "n")
lines(x, 1/x, col = 2, lwd = 2)
polygon(c(x,5), c(1/x, 5), col = "lightskyblue1")
contour(x, y, outer(x^2, y^2, FUN = "+"), add = TRUE)
# Feasible initial points V
V <- matrix(c(1, 2, 3, 3.5, .8, 4),3,2)
for(i in 1:nrow(V)){
    z <- Pf(V[i,])
    points(V[i,1], V[i,2], col = "blue", pch = 16)
    points(z[1],z[2], col = "red", pch = 16)
    arrows(V[i,1], V[i,2], z[1], z[2], length = 0.1)
}
# Infeasible initial points V
Pf <- function(v){
    g <- function(x) (x - v[1])^2 + (1/x - v[2])^2
    z <- optimize(g, c(.2,5))$min
    c(z, 1/z)
}
V <- matrix(c(.5, .4, .7, 1.7), 2, 2)
for(i in 1:nrow(V)){
    z <- Pf(V[i,])
    points(V[i,1], V[i,2], col = "blue", pch = 16)
    points(z[1],z[2], col = "red", pch = 16)
    arrows(V[i,1], V[i,2], z[1], z[2], length = 0.1)
}
dev.off()

