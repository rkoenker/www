#First shot at new monte-carlo for HK tests
#Initialize:

	n <- 100
	R <- 1000
	alphas <- c(1,.95,.9,.85)
	ns <- c(100,200)
	S <- array(0,c(length(ns),length(alphas),R,3,3))
	
#Loop: 
	for(k in 1:length(ns)){
	  n <- ns[k]
	  for(j in 1:length(alphas)){
	      alpha <- alphas[j]
	      for(i in 1:R){
		y <- filter(rt(n,4),alpha,method="recursive")
		S[k,j,i,,] <- uroot.test(y)
		}
	      }
	  }

#NB.  The array S contains a 3 by 3 matrix of test results for each cell
#
#	row one is the results of the unmodified test for the 3 score fns
#	row two is the delta^2 to use to compute the Hansen critical value
#	row three is the modified HK statistic.
