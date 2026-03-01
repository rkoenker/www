subroutine drqssbr(m,n,m5,n2,a,b,t,toler,x,e,s,wa,wb,nsol,sol,lsol,lam,swa,big,eps,dim,k0)
integer i,j,k,kl,kount,kr,l,lsol,m,m1,m2,m3,m4,m5
integer n,n1,n2,nsol,out,s(m),nobs,dim,k0
logical stage,test,init,itend,ilend,ilfix
real a(m,n),b(m),t,sol(n+6,nsol),x(n),e(m),toler,lam,t0,t1,tnt,l0,l1,lnt
real swa(m5,n2),eps,big
double precision a1,aux,b1,d,dif,pivot,smax
double precision min,max,zero,half,one,two,three,four,five
double precision wa(m5,n2),wb(m),sum,penal,fidel
data zero/0.d0/
data half/0.5d0/
data one/1.d0/
data two/2.d0/
data three/3.d0/
data four/4.d0/
data five/5.d0/
#
#  initialization
#
m1 = m+1
n1 = n+1
nobs=m-n+2
m2 = m+2
m3 = m+3
m4 = m+4
m5 = m+5
wa(m2,n1) = zero
if (n2!=n+2)
	wa(m2,n1) = three
if (m<=zero||n<=zero)
	wa(m2,n1) = four
itend = .true.
ilend = .true.
ilfix = .true.
tnt = t
if (t<zero||t>one) {
	t0 = one/dfloat(m)-toler
	t1 = one-toler  
	t = t0	
	itend = .false.
	}
if (lam<zero) {
	l0 = toler
	l1 =big-toler
	lam = l1**.25
	ilend = .false.
	ilfix = .false.
	}
if(!itend&!ilend)
	wa(m2,n1)=five #don't allow both t and lam to vary
if (wa(m2,n1)<=two) {
	do i = 1,m {
		wb(i) = b(i)
		do j = 1,n
			wa(i,j) = a(i,j) #note no lam here
		}
	dif = zero
	init = .false.
	lsol = 1
	kount = 0
	do j = 1,n {
		wa(m5,j) = j
		x(j) = zero
		}
	do i = 1,m {
		wa(i,n2) = n+i
		wa(i,n1) = wb(i)
		if (wb(i)<zero)
			do j = 1,n2
				wa(i,j) = -wa(i,j)
		e(i) = zero
		}
	do j = 1,n {
		wa(m2,j) = zero
		wa(m3,j) = zero
		wa(m4,j) = zero
		do i = 1,nobs {
			aux = dsign(one,wa(m5,j))*wa(i,j)
			wa(m2,j) = wa(m2,j)+aux*(one-dsign(one,wa(i,n2)))
			wa(m4,j) = wa(m4,j)+aux*dsign(one,wa(i,n2))
			}
		do i = nobs+1,m {
			aux = dsign(one,wa(m5,j))*wa(i,j)
			wa(m3,j) = wa(m3,j)+aux
			}
		wa(m4,j) = two*wa(m4,j)
		}
	repeat {
#
# compute new marginal costs
#
# lam--new lambda
#
		do j = 1,n{
			wa(m1,j) = wa(m2,j)+wa(m3,j)*lam+wa(m4,j)*t
			}
		if (!init) {
#
# stage 1
#
# determine the vector to enter the basis
#
			stage = .true.
			kr = 1
			kl = 1
			go to 30
			}
		repeat {
#
# stage 2
#
			stage = .false.
			repeat {
#
# determine the vector to enter the basis
#
				max = -big
				do j = kr,n {
					d = wa(m1,j)
					if (d<zero) {
						if(abs(wa(m5,j))>n+nobs){
							if (d>(-two*lam))
								next 1
							d = -d-two*lam
							}
						else {
							if (d>(-two))
								next 1
							d = -d-two
							}
						}
					if (d>max) {
						max = d
						in = j
						}
					}
				if (max<=toler)
					break 2
				if (wa(m1,in)<=zero) {
					do i = 1,m5
						wa(i,in) = -wa(i,in)
					if(abs(wa(m5,in))<=(n+nobs)){
						wa(m1,in) = wa(m1,in)-two
						wa(m2,in) = wa(m2,in)-two
						}
					else{
						wa(m1,in) = wa(m1,in)-two*lam
						wa(m3,in) = wa(m3,in)-two
						}
					}
				repeat {
#
# determine the vector to leave the basis
#
					k = 0
					do i = kl,m {
						d = wa(i,in)
						if (d>toler) {
							k = k+1
							wb(k) = wa(i,n1)/d
							s(k) = i
							test = .true.
							}
						}
					repeat {
						if (k<=0)
							test = .false.
						else {
							min = big
							do i = 1,k
								if (wb(i)<min) {
									j = i
									min = wb(i)
									out = s(i)
									}
							wb(j) = wb(k)
							s(j) = s(k)
							k = k-1
							}
#
# check for linear dependence in stage 1
#
						if (!test&&stage)
							break 1
						if (!test)
							break 5
						if(abs(wa(out,n2))>n+nobs)
							pivot = wa(out,in)*lam
						else
							pivot = wa(out,in)
#
# update wa(m1,j) and wa(m3,j) if vector leaving the basis is from the lower
# portion of the design matrix X; otherwise, update wa(m1,j) and wa(m2,j)
#
						if (wa(m1,in)-pivot-pivot<=toler)
							go to 10
						if(abs(wa(out,n2))>(n+nobs)){
							do j = kr,n1 {
								d = wa(out,j)
								wa(m1,j) = wa(m1,j)-d*lam-d*lam
								wa(m3,j) = wa(m3,j)-d-d
								wa(out,j) = -d
								}
							}
						else {
							do j = kr,n1 {
								d = wa(out,j)
								wa(m1,j) = wa(m1,j)-d-d
								wa(m2,j) = wa(m2,j)-d-d
								wa(out,j) = -d
								}
							}
						wa(out,n2) = -wa(out,n2)
						}
					do i = 1,m5 {
						d = wa(i,kr)
						wa(i,kr) = wa(i,in)
						wa(i,in) = d
						}
					kr = kr+1
					go to 20
#
# pivot on wa(out,in)
#
					10  do j = kr,n
						if (j!=in){
							wa(out,j) = wa(out,j)/wa(out,in)
							}
					wa(out,n1) = wa(out,n1)/wa(out,in)
#
# note that wa(m1,j) through wa(m4,j) are all updated
#
					do i = 1,m4
						if (i!=out) {
							d = wa(i,in)
							do j = kr,n
								if (j!=in)
									wa(i,j) = wa(i,j)-d*wa(out,j)
								wa(i,n1) = wa(i,n1)-d*wa(out,n1)
							}
					do i = 1,m4
						if (i!=out)
							wa(i,in) = -wa(i,in)/wa(out,in)
					wa(out,in) = one/wa(out,in)
					d = wa(out,n2)
					wa(out,n2) = wa(m5,in)
					wa(m5,in) = d
					kount = kount+1
					if (!stage)
						break 1
#
# interchange rows in stage 1
#
					kl = kl+1
					do j = kr,n2 {
						d = wa(out,j)
						wa(out,j) = wa(kount,j)
						wa(kount,j) = d
						}
					20  if (kount+kr==n1)
						break 2
					30  max = -one
					do j = kr,n
						if (abs(wa(m5,j))<=n) {
							d = abs(wa(m1,j))
							if (d>max) {
								max = d
								in = j
								}
							}
					if (wa(m1,in)<zero)
						do i = 1,m5
							wa(i,in) = -wa(i,in)
					}
				}
			}
		wa(m2,n1) = zero
		if (kr==1) {
			do j = 1,n {
				d = abs(wa(m1,j))
				if (d<=toler||two-d<=toler)
					go to 50
				}
			wa(m2,n1) = one
			}
		50 sol(3,lsol) = wa(m2,n1)
		do i = 1,kl-1 {
			k = wa(i,n2)*dsign(one,wa(i,n2))
			x(k) = wa(i,n1)*dsign(one,wa(i,n2))
			}
		if (itend & ilend)
			go to 40
#
# compute dimensionality, fidelity and roughness of the fit
#
		dim = 0
		penal = zero	
		fidel = zero
		do i = kl,m {
			k = wa(i,n2)*dsign(one,wa(i,n2))
			d = wa(i,n1)*dsign(one,wa(i,n2))
			k = k-n
			if(k>nobs)
				penal = penal+d*dsign(one,d)
			else if (k>0)
				fidel = fidel+d*dsign(one,d)*(one+dsign(one,d)*(two*tnt-one))
			}
		do i = 1,n{
			if (abs(wa(m5,i))-n<=nobs)
				dim = dim + 1
			}
		sol(4,lsol) = fidel
		sol(5,lsol) = penal
		sol(6,lsol) = dim
		do i = 1,n
			sol(i+6,lsol) = x(i)
		if (ilend){
			sol(1,lsol) = smax
			sol(2,lsol) = lam
			}
		if (itend){
			sol(1,lsol) = t
			sol(2,lsol) = smax
			}
		lsol = lsol+1
		init = .true.
		do i = 1,m
			s(i) = zero
		do j = 1,n
			x(j) = zero
#
#  compute next t or next lam
#
		if(ilfix){
			smax = two
			do j = 1,n {
				if(abs(wa(m5,j))>n+nobs){
					b1 = wa(m4,j)/lam
					a1 = -(one+(wa(m2,j)+wa(m3,j)*lam+lam)/lam)/b1
					b1 = (one-(wa(m2,j)+wa(m3,j)*lam+lam)/lam)/b1
					if (a1>=t)
						if (a1<smax){
							smax = a1
							dif = (b1-a1)/two
						}
					if (b1>t)
						if (b1<smax) {
							smax = b1
							dif = (b1-a1)/two
							}
					}
				else{
					b1 = wa(m4,j)
					a1 = (-two-(wa(m2,j)+wa(m3,j)*lam))/b1
					b1 = (-wa(m2,j)-wa(m3,j)*lam)/b1
					if (a1>=t)
						if (a1<smax) {
							smax = a1
							dif = (b1-a1)/two
							}
					if (b1>t)
						if (b1<smax) {
							smax = b1
							dif = (b1-a1)/two
							}
					}
				}
			tnt = smax
			t = smax+toler*(one+abs(dif))*1.d4  
			if (t>=t1+toler)  
				itend = .true.
			if (itend)
				t = t1
			}
		else {
			smax = l0
			do j = 1,n{
				if(abs(wa(m5,j))>(n+nobs)){
					b1 = wa(m3,j)
						a1 = -(wa(m2,j)+wa(m4,j)*t)/(two+b1)
						b1 = -(wa(m2,j)+wa(m4,j)*t)/b1
						}
				else{
					b1 = wa(m3,j)
						a1 = (-two-(wa(m2,j)+wa(m4,j)*t))/b1
						b1 = (-wa(m2,j)-wa(m4,j)*t)/b1
					}
				if (a1<=lam)
					if (a1>smax){
						smax = a1
						if (abs(a1)>big**.25|abs(b1)>big**.25) #avoid abnormal a1 and b1
							dif = one
						else
							dif = (b1-a1)/two
						}
				if (b1<lam)
					if (b1>smax){ 
						smax = b1
						if (abs(a1)>big**.25|abs(b1)>big**.25) #avoid abnormal a1 and b1
							dif = one
						else
							dif = (b1-a1)/two
						}
				}
			lnt = smax-toler*(one+abs(dif))*1.d4  
			if (lnt<l0|dim>=k0)
				ilend = .true.
			lam = lnt
			if (ilend){
				lam = smax+eps
				go to 40 
				}
			}
		}
	wa(m2,n1) = two
	sol(3,lsol) = wa(m2,n1)
	go to 70
#
# tidy up the sol matrix
#
	40  if (lsol>2) {
		if (ilfix){
			sol(1,lsol) = one
			sol(2,lsol) = lam
			sol(4,lsol) = zero
			sol(5,lsol) = zero
			sol(6,lsol) = two
			sol(1,1) = zero
			sol(2,1) = lam
			sol(4,1) = zero
			sol(5,1) = zero
			sol(6,1) = two
			do i=1,n
				sol(6+i,lsol) = sol(6+i,lsol-1)
			}
		else{
			sol(1,lsol) = t
			sol(2,1) = l1
			sol(2,lsol) = lam
			sol(3,lsol) = sol(3,lsol-1)
			sol(4,lsol) = sol(4,lsol-1)
			sol(5,lsol) = sol(5,lsol-1)
			sol(6,lsol) = sol(6,lsol-1)
			do i=1,n
				sol(6+i,lsol) = sol(6+i,lsol-1)
			}
		}
	l = kl-1
	do i = 1,l
		if (wa(i,n1)<zero)
			do j = kr,n2
				wa(i,j) = -wa(i,j)
	70 sum = zero
#
#compute the residuals and the norm for single t and lam
#
	do i = kl,m {
		k = wa(i,n2)*dsign(one,wa(i,n2))
		d = wa(i,n1)*dsign(one,wa(i,n2))
		k = k-n
		e(k) = d
		if(k>nobs)
			sum = sum+d*dsign(1.0,d)*lam
		else if (k>0)
			sum = sum+d*dsign(one,d)*(one+dsign(one,d)*(two*t-one))
		}
	dim = 0
	do i = 1,n
		if (abs(wa(m5,i))-n<=nobs)
			dim = dim + 1
	wa(m2,n2) = kount
	wa(m1,n2) = n1-kr
	wa(m1,n1) = sum
	}
do i=1,m5
	do j=1,n2
		swa(i,j)=wa(i,j)
return
end
