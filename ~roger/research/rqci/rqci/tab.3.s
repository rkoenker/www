#analysis phase -- Table 3 -- see monte.s for code to compute aa.i.i arrays
#alist_paste("a.",rep(c(1,3,8),3),".",rep(c(1,3,8),rep(3,3)),sep="")
coverage_matrix(0,5,3)
cilength_matrix(0,5,3)
dfx_c(1,3,8)
dfy_c(1,3,8)
methods_c("rank-inverse","sparsity-HS","PWY","Heqf-BS","XY-BS")
for(i in 1:3){
	Sum_sum.A(get(paste("aa.",dfx[i],sep="")))
	coverage[,i]_1-Sum$cov
	cilength[,i]_Sum$cil
	}
table.3_cbind(coverage, cilength)
dimnames(table.3)_list(methods,rep(paste("df=",dfy),2))
latex.table(table.3,dec=3,cgroup=c("coverage","length"),n.cgroup=c(3,3),
	rowlabel="", caption="Confidence Interval Performance -- 
	Heteroskastic Errors")
