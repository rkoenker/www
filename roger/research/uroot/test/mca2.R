#Now make latex tables 

	alphas <- c(1.00,0.95,0.90,0.85)
	ns <- c(100,200)

D <- array(c(Dtab.normal,Dtab.t2,Dtab.t3,Dtab.t4),dim=c(2,4,2,3,4))
DM <- D[,,2,,]
DU <- D[,,1,,]
anames <- paste("$\\alpha =$",alphas)
scores <- c("Wilcoxon", "Normal", "Sign")
dfs <- c("Normal", "Student on 2", "Student on 3", "Student on 4")
DMtab <- NULL
DUtab <- NULL
for(i in 1:length(dfs)){
	DMtab <- rbind(DMtab,cbind(DM[1,,,i],DM[2,,,i]))
	DUtab <- rbind(DUtab,cbind(DU[1,,,i],DU[2,,,i]))
	}
dimnames(DMtab) <- list(rep(anames,length(dfs)),rep(scores,length(ns)))
dimnames(DUtab) <- list(rep(anames,length(dfs)),rep(scores,length(ns)))
latex.table(DMtab,file = "HKm",cgroup = paste("$n =$",ns),n.cgroup = c(3,3),
	rgroup = dfs, n.rgroup = c(4,4,4,4),caption = "Modified Test")
latex.table(DUtab,file = "HKu",cgroup = paste("$n =$",ns),n.cgroup = c(3,3),
	rgroup = dfs, n.rgroup = c(4,4,4,4),caption = "Unmodified Test")
