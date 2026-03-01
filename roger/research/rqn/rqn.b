# To unbundle, sh this file
echo fnc.r 1>&2
cat >fnc.r <<'End of fnc.r'
#Outer wrapper for the new rqn function--calls a frisch-newton LP solver
#Does preprocessing using the functions globit and checkit to reduce initial n
subroutine rqm(n2,p,a,y,rhs,d,wn,wp,beta,eps,tau,s,aa,hist)
integer n,p,s(n2),hist(3,32),kit,nit,mit,m,mm,n2,maxnit,maxmit,mlim
double precision a(p,n2),y(n2),rhs(p),d(n2),wn(1),wp(p,p+3),aa(p,p)
double precision beta,eps,tau,omega,sparsity
#real ut,time,udt

data zero/0.0d0/
data one/1.0d0/
data two/2.0d0/

maxmit=8
maxnit=4
n=n2-2
m=2*nint(n**(2./3.))
mlim = 5*m
mit=0
kit=0

#outer loop on the initial sample size
while(mit<maxmit){
	mit=mit+1
	kit=kit+1
	nit=1
	call rqfn(m,p,a,y,rhs,d,beta,eps,tau,wn,wp,aa,hist(1,kit))
	omega=two*sqrt(tau*(one-tau))*sparsity(m,wn,tau)
	call dscal(p*p,omega,aa,1)
	call globit(p,n2,aa,a,y,s,wp,wp(1,2),wp(1,3),wp(1,4),mm)
	#inner loop to check for optimality and fixup a few bad globbed apples
	while(nit<maxnit){
		m=mm
		nit=nit+1
		kit=kit+1
		call rqfn(mm,p,a,y,rhs,d,beta,eps,tau,wn,wp,aa,hist(1,kit))
		call checkit(n2,p,a,y,s,wp,mm,mlim)
		#try to fixup or, if too bad, draw another initial sample
		if(mm>m){if(mm>n){m=2*m;break}}
		else return
		}
	}
return
end
#timing function
real function udt(t)
real s,dtime,udt,u(2)
s=dtime(u)
udt=u(1)
return
end

#This is my second attempt to code the primal-dual log barrier form of the
#interior point LP solver of Lustig, Marsten and Shanno ORSA J Opt 1992.
#It is a projected Newton primal-dual logarithmic barrier method which uses
#the predictor-corrector approach of Mehrotra for the mu steps.
#For the sake of brevity we will call it a Frisch-Newton algorithm.
#The primary difference between this code and the previous version fna.r is
#that we assume a feasible starting point so p,d,b feasibility gaps = 0.
#Note also that all variables are assumed to have upper bounds of one.
#Problem:
#	min c'x s.t. Ax=b, 0<=x<=1
#
#The linear system we are trying to solve at each interation is:
#		 Adx = 0
#	     dx + ds = 0
#	A'dy -dw +dz = 0
#	   Xdz + Zdx = me - XZe - DXDZe
#	   Sdw + Wds = me - SWe - DSDWe 
#But the algorithm proceeds in two steps the first of which is to solve:
#		 Adx = 0
#	     dx + ds = 0
#	A'dy -dw +dz = 0
#	   Xdz + Zdx =  - XZe 
#	   Sdw + Wds =  - SWe
#and then to make some refinement of mu and modify the implied Newton step.
#Denote dx,dy,dw,ds,dz as the steps for the respective variables x,y,w,s,z, and
#by the corresponding upper case letters are the diagonal matrices respectively.
#
#To illustrate the use of the function we include a calling routine to
#compute solutions to the linear quantile regression estimation problem.
#See the associated S function rqfn for further details on the calling sequence.
subroutine rqfn(n,p,a,y,rhs,d,beta,eps,tau,wn,wp,aa,nit)
integer n,p,nit(3)
double precision a(p,n),y(n),rhs(p),d(n),beta,eps,wn(n,10),wp(p,p+3),aa(p,p)
double precision mone,one,tau,ddot
data one/1.0d0/
data mone/-1.0d0/
do i=1,n{
	d(i)=one
	wn(i,1)=one-tau
	}
do i=1,p
	rhs(i)=(one-tau)*ddot(n,d,1,a(i,1),p)
call fna(n,p,a,y,rhs,d,beta,eps,wn(1,1),wn(1,2),
	wp(1,1),wn(1,3),wn(1,4),wn(1,5), wn(1,6),
	wp(1,2),wn(1,7),wn(1,8),wn(1,9),wn(1,10),wp(1,3), wp(1,4),aa,nit)
return
end
subroutine fna(n,p,a,c,b,d,beta,eps,x,s,y,z,w,
		dx,ds,dy,dz,dw,dsdw,dxdz,rhs,ada,aa,nit)

integer n,p,pp,i,info,nit(3)
double precision a(p,n),c(n),b(p)
double precision zero,one,mone,big,ddot,dmax1,dmin1,dasum
double precision deltap,deltad,beta,eps,cx,by,uw,uz,mu,mua,acomp,rdg,g
double precision x(n),s(n),y(p),z(n),w(n),d(n),rhs(p),ada(p,p),aa(p,p)
double precision dx(n),ds(n),dy(p),dz(n),dw(n),dxdz(n),dsdw(n)

data zero /0.0d0/
data half /0.5d0/
data one  /1.0d0/
data mone  /-1.0d0/
data big /1.0d+20/

#Initialization:  We try to follow the notation of LMS
#On input we require:
#
#	c=n-vector of marginal costs (-y in the rq problem)
#	a=p by n matrix of linear constraints (x' in rq)
#	b=p-vector of rhs ((1-tau)x'e in rq)
#	beta=barrier parameter, LMS recommend .99995
#	eps =convergence tolerance, LMS recommend 10d-8
#	
nit(1)=0
nit(2)=0
nit(3)=n
pp=p*p
#Start at the OLS estimate for the parameters
call dgemv('N',p,n,one,a,p,c,1,zero,y,1)
call stepy(n,p,a,d,y,aa)
#Save sqrt of aa' for future use for confidence band
do i=1,p{
	do j=1,p
		ada(i,j)=zero
        ada(i,i)=one
	}
call dtrtrs('U','T','N',p,p,aa,p,ada,p,info)
call dcopy(pp,ada,1,aa,1)
#put current residual vector in s (temporarily)
call dcopy(n,c,1,s,1)
call dgemv('T',p,n,mone,a,p,y,1,one,s,1)
#Initialize remaining variables
#N.B. x must be initialized on input: for rq as (one-tau) in call coordinates
do i=1,n{
	d(i)=one
	z(i)=dmax1(s(i),zero)
	w(i)=dmax1(-s(i),zero)
	s(i)=one-x(i)
	}
cx = ddot(n,c,1,x,1)
by = ddot(p,b,1,y,1)
uw = dasum(n,w,1)
uz = dasum(n,z,1)
#rdg =  (cx - by + uw)/(one + dabs( by - uw))
#rdg =  (cx - by + uw)/(one + uz + uw)
rdg =  (cx - by + uw)
while(rdg > eps) {
	nit(1)=nit(1)+1
	do i =1,n{
		d(i) = one/(z(i)/x(i) + w(i)/s(i))
		ds(i)=z(i)-w(i)
		dx(i)=d(i)*ds(i)
		}
	call dgemv('N',p,n,one,a,p,dx,1,zero,dy,1)#rhs
	call dcopy(p,dy,1,rhs,1)#save rhs
	call stepy(n,p,a,d,dy,ada)
	call dgemv('T',p,n,one,a,p,dy,1,mone,ds,1)
	deltap=big
	deltad=big
	do i=1,n{
		dx(i)=d(i)*ds(i)
		ds(i)=-dx(i)
		dz(i)=-z(i)*(dx(i)/x(i) + one)
		dw(i)=w(i)*(dx(i)/s(i) - one)
		dxdz(i)=dx(i)*dz(i)
		dsdw(i)=ds(i)*dw(i)
		if(dx(i)<0)deltap=dmin1(deltap,-x(i)/dx(i))
		if(ds(i)<0)deltap=dmin1(deltap,-s(i)/ds(i))
		if(dz(i)<0)deltad=dmin1(deltad,-z(i)/dz(i))
		if(dw(i)<0)deltad=dmin1(deltad,-w(i)/dw(i))
		}
	deltap=dmin1(beta*deltap,one)
	deltad=dmin1(beta*deltad,one)
	if(deltap*deltad<one){
		nit(2)=nit(2)+1
		acomp=ddot(n,x,1,z,1)+ddot(n,s,1,w,1)
		g=acomp+deltap*ddot(n,dx,1,z,1)+
			deltad*ddot(n,dz,1,x,1)+ 
			deltap*deltad*ddot(n,dz,1,dx,1)+
			deltap*ddot(n,ds,1,w,1)+
			deltad*ddot(n,dw,1,s,1)+ 
			deltap*deltad*ddot(n,ds,1,dw,1)
		mu=acomp/dfloat(2*n)
		mua=g/dfloat(2*n)
		mu=mu*(mua/mu)**3
		#if(acomp>1) mu=(g/dfloat(n))*(g/acomp)**2
		#else mu=acomp/(dfloat(n)**2)
		do i=1,n{
			dz(i)=d(i)*(mu*(1/s(i)-1/x(i))+
				dx(i)*dz(i)/x(i)-ds(i)*dw(i)/s(i))
			}
		call dswap(p,rhs,1,dy,1)
		call dgemv('N',p,n,one,a,p,dz,1,one,dy,1)#new rhs
		call dpotrs('U',p,1,ada,p,dy,p,info)
		call daxpy(p,mone,dy,1,rhs,1)#rhs=ddy
		call dgemv('T',p,n,one,a,p,rhs,1,zero,dw,1)#dw=A'ddy
		deltap=big
		deltad=big
		do i=1,n{
			dx(i)=dx(i)-dz(i)-d(i)*dw(i)
			ds(i)=-dx(i)
			dz(i)=mu/x(i) - z(i)*dx(i)/x(i) - z(i) - dxdz(i)/x(i)
			dw(i)=mu/s(i) - w(i)*ds(i)/s(i) - w(i) - dsdw(i)/s(i)
			if(dx(i)<0)deltap=dmin1(deltap,-x(i)/dx(i))
			else deltap=dmin1(deltap,-s(i)/ds(i))
			if(dz(i)<0)deltad=dmin1(deltad,-z(i)/dz(i))
			if(dw(i)<0)deltad=dmin1(deltad,-w(i)/dw(i))
			}
		deltap=dmin1(beta*deltap,one)
		deltad=dmin1(beta*deltad,one)
		}
	call daxpy(n,deltap,dx,1,x,1)
	call daxpy(n,deltap,ds,1,s,1)
	call daxpy(p,deltad,dy,1,y,1)
	call daxpy(n,deltad,dz,1,z,1)
	call daxpy(n,deltad,dw,1,w,1)
	cx=ddot(n,c,1,x,1)
	by=ddot(p,b,1,y,1)
	uw = dasum(n,w,1)
	uz = dasum(n,z,1)
	#rdg=(cx-by+uw)/(one+dabs(by-uw))
	#rdg=(cx-by+uw)/(one+uz+uw)
	rdg=(cx-by+uw)
	}
#return residuals in the vector x
call daxpy(n,mone,w,1,z,1)
call dswap(n,z,1,x,1)
return
end
subroutine stepy(n,p,a,d,b,ada)
integer n,p,pp,i,info
double precision a(p,n),b(p),d(n),ada(p,p),zero
data zero/0.0d0/
#Solve linear system ada'x=b by Choleski -- d is diagonal
#Note that a isn't altered, and on output ada is the 
#upper triangle Choleski factor, which can be reused, eg with blas dtrtrs
pp=p*p
do j=1,p
	do k=1,p
		ada(j,k)=zero
do i=1,n
	call dsyr('U',p,d(i),a(1,i),1,ada,p)
call dposv('U',p,1,ada,p,b,p,info)
return
end
End of fnc.r
echo glob.r 1>&2
cat >glob.r <<'End of glob.r'
#This subroutine reorders the sample and makes the globs
#On input:
#	aa is a p by p upper choleski factor of x'x inverse times omega
#			i.e.  the estimated sqrt of the cov matrix of bhat
#	a  is a p by n+2 design matrix (note transpose!)
#	b  is the observed response again with n+2 elements to accomodate globs
#	bhat is fitted coefficient vector
#	
#On output
#	aa is unaltered
#	a  is reordered with first nxt columns for globbed sample
#		and remaining columns to save the  globbed observations
#	b  is similarly reordered
#	nxt is the number of observations in the globbed sample 
subroutine globit(p,n2,aa,a,b,s,bhat,glob,ghib,work,nxt)
integer n,n2,p,lxt,nxt,s(n2)
double precision aa(p,p),a(p,n2),b(n2),bhat(p),glob(p),ghib(p),work(p)
double precision one, zero,resid,band,dnrm2,ddot,blo,bhi

data zero/0.0d0/
data one/1.0d0/

n=n2-2
blo=zero
bhi=zero
nxt=1
lxt=n
do i=1,p{
	glob(i)=zero
	glob(i)=zero
	}
while(nxt <= lxt){
	resid=b(nxt)-ddot(p,a(1,nxt),1,bhat,1)
	call dgemv('N',p,p,one,aa,p,a(1,nxt),1,zero,work,1)
	band=dnrm2(p,work,1)
	if(resid < band){
		if(resid > -band){ #inside the band points
			s(nxt)=0
			nxt=nxt+1
			}
		else{ #below the band points
			call daxpy(p,one,a(1,nxt),1,glob,1)
			call dswap(p,a(1,nxt),1,a(1,lxt),1)
			blo=blo+b(nxt)
			call dswap(1,b(nxt),1,b(lxt),1)
			s(lxt)=-1
			lxt=lxt-1
			}
		}
	else{ #above the band points
		call daxpy(p,one,a(1,nxt),1,ghib,1)
		call dswap(p,a(1,nxt),1,a(1,lxt),1)
		bhi=bhi+b(nxt)
		call dswap(1,b(nxt),1,b(lxt),1)
		s(lxt)=1
		lxt=lxt-1
		}
	}
#At this point we swap-in the globs just after the "inside" points
call dswap(p,a(1,nxt),1,a(1,n+1),1)
call dswap(p,a(1,nxt+1),1,a(1,n+2),1)
call dswap(1,b(nxt),1,b(n+1),1)
call dswap(1,b(nxt+1),1,b(n+2),1)
call dcopy(p,glob,1,a(1,nxt),1)
call dcopy(p,ghib,1,a(1,nxt+1),1)
s(n+1)=s(nxt)
s(n+2)=s(nxt+1)
nxt=nxt+1
b(nxt-1)=blo
b(nxt)=bhi
return
end
#This subroutine checks to see whether the globs are ok, if not it adjusts them.
#On input:
#	a is x' (p by n)
#	y is the response
#	s is an indicator vector: 1 in ghib, -1 in glob, 0 in sample
#	bhat is current parameter estimate based on globbed sample
#	m is the number of obs in the globbed sample
#	mlim is the maximum number of allowed adjusted observations
#On output:
#	a is the reordered x' matrix
#	y is the reordered y vector
#	m is the number of observations in the adjusted sample
#		(if m=n2 then too many fixups were encountered and the globs
#		 were exiled to the end of the full data set, in this case
#		 the algorithm should start again with a new initial sample)
subroutine checkit(n2,p,a,y,s,bhat,m,mlim)
integer n2,p,s(n2),m,nbad,mlim
double precision a(p,n2),y(n2),bhat(p)
double precision mone,resid,ddot
data mone/-1.0d0/
nbad=0
nxt=m+1
do i=m+1,n2{
	if(nxt>mlim)break
	resid=y(i)-ddot(p,a(1,i),1,bhat,1)
	if(s(i)*resid>0)next
	else {
		if(s(i)<0){#glob problem
			call daxpy(p,mone,a(1,i),1,a(1,m-1),1)
			y(m-1)=y(m-1)-y(i)
			}
		else{#ghib problem
			call daxpy(p,mone,a(1,i),1,a(1,m),1)
			y(m)=y(m)-y(i)
			}
		call dswap(p,a(1,i),1,a(1,nxt),1)
		call dswap(1,y(i),1,y(nxt),1)
		call iswap(1,s(i),s(nxt))
		nxt=nxt+1
		}
	}
if(nxt>m+1){ #swap the globs to the end of the "good" data
	if(nxt>mlim)nxt=n2+1 #abort fixup phase swap globs to (n+1,n+2)-land
	call dswap(p,a(1,m-1),1,a(1,nxt-2),1)
	call dswap(p,a(1,m),1,a(1,nxt-1),1)
	call dswap(1,y(m-1),1,y(nxt-2),1)
	call dswap(1,y(m),1,y(nxt-1),1)
	call iswap(1,s(m-1),s(nxt-2))
	call iswap(1,s(m),s(nxt-1))
	m=nxt-1
	}
return
end
subroutine iswap(n,a,b)
#Integer swap routine to mimic lapack dswap
integer n, a(n),b(n),aa
for(i=1;i<=n;i=i+1){
	aa=a(i)
	a(i)=b(i)
	b(i)=aa
	}
return
end
End of glob.r
echo sparsity.r 1>&2
cat >sparsity.r <<'End of sparsity.r'
#This is a Siddiqui sparsity function estimate based on residuals
double precision function sparsity(n,u,tau)
integer n,nd,enuf
double precision u(n),tau,h,qhi,qlo,half

data half/0.5d0/
data enuf/600/

#bandwidth: approximate Hall-Sheather method - quadratic approx max error 1%
nd=nint((.05 + 3.65*tau - 3.65*tau**2)*(n**(2./3.)))
h=dfloat(nd)/dfloat(n)

#compute Siddiqui estimator
call kuantile(n,u,tau+h,qhi,enuf,half,half)
call kuantile(n,u,tau-h,qlo,enuf,half,half)
sparsity=(qhi-qlo)/(h+h)
return
end
#function to compute pth quantile of a sample of n observations
subroutine kuantile(n,x,p,q,mmax,cs,cd)
integer n,k,l,r,mmax
double precision  x(n),p,q,cs,cd
if(p<0 | p>1) {call dblepr("sparsity bandwidth problem: p=",30,p,1);return}
l=1
r=n
k=nint(p*n)
call select(n,x,l,r,k,mmax,cs,cd)
q=x(k)
return
end
#This is a ratfor implementation of the floyd-revest quantile algorithm--SELECT
#Reference:  CACM 1975, alg #489, p173, algol-68 version
#As originally proposed:  mmax=600, and cs=cd=.5
#Calls blas routine dswap 
subroutine select(n,x,l,r,k,mmax,cs,cd)
integer n,m,l,r,k,ll,rr,i,j,mmax
double precision  x(n),z,s,d,t,cs,cd
while(r>l){
	if(r-l>mmax){
		m=r-l+1
		i=k-l+1
		fm=dfloat(m)
		z=log(fm)
		s=cs*exp(2*z/3)
		d=cd*sqrt(z*s*(m-s)/fm)*sign(1.,i-m/2)
		ll=max(l,k-i*s/fm +d)
		rr=min(r,k+(m-i)*s/fm +d)
		call select(n,x,ll,rr,k,mmax,cs,cd)
		}
	t=x(k)
	i=l
	j=r
	call dswap(1,x(l),1,x(k),1)
	if(x(r)>t)call dswap(1,x(r),1,x(l),1)
	while(i<j){
		call dswap(1,x(i),1,x(j),1)
		i=i+1
		j=j-1
		while(x(i)<t)i=i+1
		while(x(j)>t)j=j-1
		}
	if(x(l)==t)
		call dswap(1,x(l),1,x(j),1)
	else{
		j=j+1
		call dswap(1,x(j),1,x(r),1)
		}
	if(j<=k)l=j+1
	if(k<=j)r=j-1
	}
return
end
End of sparsity.r
echo makefile 1>&2
cat >makefile <<'End of makefile'
#This is a Makefile for the new frisch-newton RQ routine
#These flags are intended to optimize for ysidro, to speed compile drop them
CFLAGS = -c -xarch=v8 -xchip=ultra -O4
#These flags are intended to optimize for ragnar, to speed compile drop them
#CFLAGS = -c -xarch=v8 -xchip=super2 -O4
#CFLAGS = -c 
LFLAGS = -r -dn /usr/local/SUNWspro/SC4.0/lib/v8/libsunperf.a

fn.o: fnc.o glob.o sparsity.o
	ld fnc.o glob.o sparsity.o $(LFLAGS) -o fn.o
fnc.o: fnc.r
	f77 $(CFLAGS) fnc.r
glob.o: glob.r
	f77 $(CFLAGS) glob.r
sparsity.o: sparsity.r
	f77 $(CFLAGS) sparsity.r
fb.o: fnb.o
	ld fnb.o $(LFLAGS) -o fb.o
fnb.o: fnb.r
	f77 $(CFLAGS) fnb.r
fa.o: fna.o
	ld fna.o $(LFLAGS) -o fa.o
fna.o: fna.r
	f77 $(CFLAGS) fna.r
clean:
	rm fnc.o fn.o

End of makefile
