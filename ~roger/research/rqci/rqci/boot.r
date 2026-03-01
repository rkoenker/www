#rq for multiple y's
subroutine rqs(m,n,k,m5,n2,a,b,t,toler,ift,x,e,s,
		wa,wb,nsol,ndsol,sol,dsol,lsol,h)
real b(m,k),sol(n2,nsol),a(m,n),x(n,k),wa(m5,n2),wb(m), e(m),dsol(m,nsol)
integer h(n,nsol),s(m)
do i=1,k 
	call rq(m,n,m5,n2,a,b(1,i),t,toler,ift,x(1,i),e,s,
		wa,wb,nsol,ndsol,sol,dsol,lsol,h)
return
end
#parzen, wei and ying's bootstrap
subroutine pwy(m,n,k,m5,n2,a,c,b,t,toler,ift,x,e,s,
		wa,wb,nsol,ndsol,sol,dsol,lsol,h)
real b(m),sol(n2,nsol),a(k,n),x(n,k),wa(m5,n2),wb(m),e(m),dsol(m,nsol),c(m,n)
integer h(n,nsol),s(m)
do i=1,k{
	call scopy(n,a(i,1),k,c(m,1),m)
	call rq(m,n,m5,n2,c,b,t,toler,ift,x(1,i),e,s,
		wa,wb,nsol,ndsol,sol,dsol,lsol,h)
	}
return
end
#ratfor outer loop for xy-pairs rq bootstrap
subroutine xys(m,n,k,m5,n2,a,b,t,toler,ift,x,e,s,
	wa,wb,nsol,ndsol,sol,dsol,lsol,h,aa,bb,ss)
real b(m),sol(n2,nsol),a(m,n),x(n,k),wa(m5,n2),wb(m),e(m),dsol(m,nsol)
real aa(m,n),bb(m)
integer h(n,nsol),ss(m,k),s(m)
do i=1,k {
	do ii=1,m{
		bb(ii)=b(ss(ii,i))
		do jj=1,n{
			aa(ii,jj)=a(ss(ii,i),jj)
			}
		}
	call rq(m,n,m5,n2,aa,bb,t,toler,ift,x(1,i),e,s,
			wa,wb,nsol,ndsol,sol,dsol,lsol,h)
	}
return
end
#does a matrix multiply to make Y matrix for heqf bootstrap
subroutine heqfy(n,p,r,x,b,y)
integer n,p,r
real x(n,p),b(p,n,r),y(n,r)
do i=1,r{
	do j=1,n{
		y(j,i)=sdot(p,x(j,1),n,b(1,j,i),1)
		}
	}
return
end
