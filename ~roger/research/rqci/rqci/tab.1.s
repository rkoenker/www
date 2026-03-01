#analysis phase
#alist_paste("a.",rep(c(1,3,8),3),".",rep(c(1,3,8),rep(3,3)),sep="")
coverage_array(0,c(3,3,3))
cilength_array(0,c(3,3,3))
dfx_c(1,3,8)
dfy_c(1,3,8)
methods_c("rank-inverse","sparsity-HS","sparsity-BS")
for(i in 1:3){
	for(j in 1:3){
		Sum_sum.A(get(paste("a.",dfx[i],".",dfy[j],sep="")))
		coverage[,i,j]_1-Sum$cov
		cilength[,i,j]_Sum$cil
		}
	}
table.1_cbind(matrix(coverage,9,3), matrix(cilength,9,3))
dimnames(table.1)_list(rep(methods,3),rep(paste("dfy=",dfy),2))
latex.table(table.1,dec=3,cgroup=c("coverage","length"),n.cgroup=c(3,3),
	rowlabel="",rgroup=paste("dfx=",dfx),n.rgroup=c(3,3,3),label="",
	caption="Confidence Interval Performance -- IID Errors")
