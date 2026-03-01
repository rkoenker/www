# Analysis of the q2 simulation
require(Hmisc)
load("q2.Rda")
A <- array(unlist(AA),c(length(ns),length(betas),R,length(lambdas)))
tab <- apply(A > cvs, c(1,2,4),mean)
tab <- aperm(tab,c(1,3,2))
tab <- matrix(tab,length(ns))
cgrp <- paste("$\\beta_0 = $", 0:2)
rnames <- paste("n = ", ns)
cnames <- rep(paste("$\\lambda = $", lambdas),3)
dimnames(tab) <- list(rnames,cnames)
cap <- paste("Rejection frequencies for the Hotelling likelihood ratio test
for a simple Box-Cox example. Tests are nominal level $\\alpha = 0.05$.   
Local alternatives are employed of the form:  $\\beta_n = \\beta_0 / \\sqrt{n}.$")
latex(tab, file = "q2.tex", rowlabel = "", cgroup = cgrp, 
	where = "!htbp", caption = cap)

