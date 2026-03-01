//*********************************************************************
// RQ Ox  version 1.0 August 1999
// This is a group of Ox programs to solve unconstrained quantile regression 
// problems using a preprocessing approach described in Portnoy and Koenker
// (1997) and  interior point methods. The Ox translation of the interior point
// algorithm is by Daniel Morillo and Roger Koenker of Splus and ratfor code by R.
// Koenker. The implementation of the parametric programming and inference procedures
// is by Daniel Morillo.
// See Portnoy and Koenker (1997) for a detailed description of the algorithm.
// The program now includes two functions that solve parametric programming
// problems once the regression is computed. This allows the computation of
// coefficients for a range of tau as well as inference on the obtained
// coefficients
//
// The programs in this release comprise the following:
//
//	lp_fnm: Frisch-Newton-Mehrotra algorithm for solving bounded
//		variables linear programs.
//	rq_fnm: Calling routine for using lp_fnm to solve rq problems.
//	sparsity: A simple function to compute Siddiqui estimate of
//		the sparsity function (reciprocal of the density). 
//	RQ: This function implements the algorithm described in
//		Portnoy and Koenker (1997).  It calls rq_fnm after
//		some preprocessing to reduce the effective dimension of
//		the problem.
//	rq_alltau: A function that solves the problem for all possible taus
//	pprog_rhs: This function does righ-hand side parametric programming
//	rq_confidence: This function first solves the problem by calling RQ
//		and then obtains confidence interval on the coefficients by
//		calling pprog_rhs
//  kerncv: This function computes a kernel-type variance covariance matrix
//		as in Powell (1996)
//	interp: This function does simple interpolation so that interpolated
//		rank-based confidence interval boundaries can be obtained
//  rq_inference: This function implements some of the inference procedures
//		in Koenker and Machado (1999)
//	rq: This is the "wrapper function" that the end user should use.
//
//
//	NOTE: The end user should be able to use just the function "rq" for most
//		purposes. All functions and global variables that are used to change
//		varios options are fully documented in the code.
//
// This version of the functions were tested in Version 2.1 of Ox. See
// http://www.nuff.ox.ac.uk/Users/Doornik/ for details on Ox.  The paper
// by Portnoy and Koenker and further details of the Splus/Fortran
// implementation of the algorithm are available from 
// http://www.econ.uiuc.edu/~roger/research/rqn/rqn.html.

//**********************************************************************

#include <oxstd.h>
#include <oxprob.h>

//*****************************************************************************
// GLOBAL VARIABLES:
//
// KCV, GLOB, RCI, ALPHA, REPORT_R1: See function rq.
//
// DUAL_SAVEMODE: See the function rq_alltau.
//
// BW_EXP, BW_CONS, KERNEL:  See the function kcv.
//
// VERBOSE: Controls output during program operation. Set to zero for no
//	intermediate output. Up to 3 for maximum verbosity.
//*****************************************************************************

decl KCV=TRUE;
decl GLOB=FALSE;
decl RCI=TRUE;
decl ALPHA=0.05;
decl REPORT_R1=TRUE;
decl REPORT_LK=TRUE;
decl DUAL_SAVEMODE=2;
decl BW_EXP=2/5;
decl BW_CONS=0;
decl KERNEL="normal";
decl VERBOSE=0;

//**********************************************************************
// Function: lp_fnm
//
// Description:
// This function is designed to solve the bounded variables linear program:
//
//   min_x {c'x | Ax = b,  0 <= x <= u}
//
// Data:
//
//   A  is a p by n matrix
//   c  is a n vector
//   b  is a p vector
//   u  is a n vector
//
// Control Parameters:
//
//   eps is a scalar termination criterion parameter
//   beta is a scalar step length parameter
//
// VERSION NOTES:
//
//  1) This is a Frisch-Newton-Mehrotra algorithm for solving bounded variables
//    linear programs.  This designation is based on the following attributions
//    Frisch(1956) who apparently first suggested log-barrier methods for LP's
//    Newton(?) who apparently first suggested Newton's method
//    Mehrotra(1992) who apparently first suggested the particular updating 
//    used here for the barrier parameter, mu.
//  2) Many others obviously had a hand in developing the ideas underlying 
//    the algorithm, among whom we might mentioned Kantorovich's
//    analysis of Newton's method, and Fiacco and McCormick's
//    development of barrier methods.  Many details of the 
//    implementation used here are based on the implementation
//    of Lustig, Marsden and Shanno(1994).
//  3) The Ox version is based on a translation of an Splus/ratfor
//    version written by R. Koenker, and described in Portnoy
//    and Koenker (1997).  The Ox translation was carried out by
//    Daniel Morillo and Roger Koenker. 
//  4) Contrary to LMS, we assume that an initial _feasible_ solution
//    is available.  This simplifies somewhat the rhs's.
//  5) The updating of mu reverts to the original Mehrotra proposal.
//  6) Convergence is based on the unnormalized duality gap.
//  7) Equation numbers below refer to LMS (1994).
//  8) All vectors are _row_ vector for annoying Ox reasons of efficiency.
//    except for dy which needs to be a column vector....aaargh!
//
// Syntax:
//
//	Input:
//		A, p by n matrix
//  	c, n row vector
//  	b, p row vector
//  	u, n row vector of upper bounds on the dual x
//  	x, address:
//			in: starting values for dual n row vector
//			out: dual n row vector
//
//	Output:
//		p column vector of primal variables at the solution
//		
//**********************************************************************
lp_fnm(const A,const c,const b,const u,x){
	decl s,w,z,y,r,q,dx,ds,dw,dz,dy,xinv,sinv,xi,dxdz,dsdw;
	decl mu,gap,delta_p,delta_d,g;
	decl rhs,L,D;
	decl it,n;
			
	decl BIG =   1e+20; //a large number
	decl MAX_IT   = 50;	//maximum number of iterations
	decl eps =   1e-05;	//convergence precision
	decl beta = 0.9995;	//step length factor

	n = columns(A);
	
	// generate initial feasible point

	if(VERBOSE>=3){
		print("\n      Starting interior point algorithm\n");
	}

	s = u-x[0];				

	ols2r(c,A,&y);			// Major Choleski's OLS

	r = c - y*A;			// residual from OLS regression
	z =  r .> 0 .? r .: 0;		// z is + part of r, w is - part of r
	w =  z - r;			// r=z-w this guarantees feasibility 

	gap = (c*x[0]'-b*y'+u*w');		// duality gap

	// loop starts	

	it  =  0;
	while( gap > eps && it < MAX_IT ){
		
		++it;
		if(VERBOSE>=3){
			print("        Doing iteration #",it," gap is ",double(gap),"\n");
		}
				
		// solve  (19) for affine step

		q = (1)./(z./x[0]+w./s);
		r = z-w;
		rhs = (q.*r)*A';
		decldl((A.*q)*A',&L,&D);// Choleski decomposition
		dy = solveldl(L,D,rhs');// backsolve for dy, NB. dy is p by 1

		//compute the other increments

		dx = q.*(dy'*A-r);
		ds = -dx;
		dz = (-z.*dx)./x[0] - z;
		dw = (-w.*ds)./s - w;

		//compute step sizes, equation (20)

		delta_p =  min((dx .< 0 .? (x[0]./(-dx)) .: BIG),
				(ds .< 0 .? (s./(-ds)) .: BIG));
		delta_d =  min((dz .< 0 .? (z./(-dz)) .: BIG),
				(dw .< 0 .? (w./(-dw)) .: BIG));

		delta_p = min(beta*delta_p,1);
		delta_d = min(beta*delta_d,1);

		//If full affine step is feasible take it, otherwise modify it.

		if(min(delta_p,delta_d)<1){

			// update mu

			mu = x[0]*z' + s*w';
			g = (x[0]+delta_p*dx)*(z+delta_d*dz)' + 
				(s+delta_p*ds)*(w+delta_d*dw)';
			mu = ((g/mu)^3)*(mu/(2*n));

			// now solve (24) for modified step

			dxdz = dx.*dz;
			dsdw = ds.*dw;
			xinv = (1)./x[0];
			sinv = (1)./s;

			xi = mu*(xinv-sinv);
			rhs = rhs+(q.*((dxdz-dsdw)-xi))*A';
			dy = solveldl(L,D,rhs');//backsolve with new rhs

			dx = q.*(dy'*A+xi-r-(dxdz-dsdw));
			ds = -dx;
			dz = mu*xinv-z-xinv.*z.*dx-dxdz;
			dw = mu*sinv-w-sinv.*w.*ds-dsdw;

			// now compute new step lengths, again as in (24)

			delta_p =  min((dx .< 0 .? (x[0]./(-dx)) .: BIG),
					(ds .< 0 .? (s./(-ds)) .: BIG));
			delta_d =  min((dz .< 0 .? (z./(-dz)) .: BIG),
					(dw .< 0 .? (w./(-dw)) .: BIG));

			delta_p = min(beta*delta_p,1);
			delta_d = min(beta*delta_d,1);
		}

		// now update all variables

		x[0] = x[0]+delta_p*dx;
		s = s+delta_p*ds;
		y = y+delta_d*dy';
		w = w+delta_d*dw;
		z = z+delta_d*dz;

		gap = (c*x[0]'-b*y'+u*w');
	}
	if(VERBOSE>=3){
		print("\n");
	}
	return(y');
}

//*********************************************************************
// Function: rq_fnm
//
// Description:
// Solves unconstrained quantile regression problems using a primal-dual
// interior point method.  Translation by Daniel Morillo and Roger Koenker of
// splus and ratfor code by R. Koenker.
//
// Syntax:
//
//	Input:
//
//		y, n column vector of n dependent observations
//		X, n by p matrix of n observations of p variables. No constant is included
//   		tau, scalar between 0 and 1 determining the quantile to estimate
//	 	a, address:
//			in: nothing
//			out: dual n row vector at the solution
//
//
// 	Output:
//   	p column vector of parameter estimates
//
//*******************************************************************
rq_fnm(const X, const y, const tau,const a){
	decl n,p,b,u;

	n = rows(X);		// sample size
	p = columns(X);		// dimensionality of beta
	u = ones(1,n);		// fixed upper bounds for dual variables
	a[0] = (1-tau)*u;	// starting values for dual variables

	b = lp_fnm(X',-y',a[0]*X,u,a); //a needs to be address
	a[0]=a[0]';
	
	return(-b);
}

//*******************************************************************
// Function: sparsity
//
// Description:
// Siddiqui estimate of the sparsity function
// A simple difference quotient estimate of the reciprocal of the
// density at a particular quantile.  This function is used to
// estimate a nuisance parameter for quantile regression inference.
//
// Syntax:
//
//	Input:
//  	r, column vector of residuals 
//  	tau, scalar determining the quantile of interest
//
//	Output:
//  	scalar, the sparsity estimate
//
//*******************************************************************
sparsity(const r, const tau){
        decl  h,qhat;
        decl n;
 
        n=rows(r);
        h=(0.05+3.65*tau-3.65*(tau^2))*(n^(2/3))/n; // Hall-Sheather bandwidth 
        qhat=quantilec(r,tau + h .* <-1 , 1>);
        return((qhat[1][0]-qhat[0][0])/(2*h)); 
}

//*********************************************************************
// Function: RQ
//
// Description:
// Solves unconstrained quantile regression problems using a primal-dual
// interior point method and the preprocessing approach described in
// Portnoy and Koenker(1997).  The Ox code is a translation by Daniel 
// Morillo and Roger Koenker of Splus and ratfor code by R. Koenker.
//
// Syntax:
//
// Input:
//   y, n column vector of independent observations
//   X, n by p matrix of n dependent observations in p variables. No constant is included
//   tau, scalar between 0 and 1 determining the quantile to estimate
//
// Output:
//   p column vector of parameter estimates
//
//*******************************************************************

RQ(X,y, const tau){
	decl n,p,m,it;
	decl not_optimal,bad_fixups;
	decl XX,yy,b,r,hi,d;
	decl su,sl,in,su_bad,sl_bad;
	decl a,S,Sinv,cs,qb,discr,b0;
	decl sortvec;

	decl MAX_BAD_FIXUPS = 5;
	decl MAX_IT = 5;

	n = rows(X);            // sample size
        p = columns(X);         // dimensionality of beta
	m = floor(max(n/20,n^(2/3)));	// preliminary sample size
	not_optimal = TRUE;
	it = 0;

	if(VERBOSE>=1){
		print("\nStarting Globbing procedure\n");
	}
	
	sortvec=rann(n,1);
	X=sortbyc(sortvec~X,0)[][1:p];
	y=sortbyc(sortvec~y,0)[][1];

	while(not_optimal && it < MAX_IT){
		++it;
		m  = 2*m;

		if(VERBOSE>=2){
			print("  Now doing iteration #",it," with ",m," observations\n");
		}

		XX = X[0:m][];
		yy = y[0:m][];
		b  = rq_fnm(XX,yy,tau,&a);	//preliminary estimation
		r  = y - X*b;
		cs = 2*sparsity(r,tau)*sqrt(tau*(1-tau));
		Sinv=XX'*XX;
		S=invertsym(Sinv);
		hi = sqrt(sumsqrr(X*choleski(S)));
		d  = cs*hi;
		su = r .>  d; 		//above the band
		sl = r .< -d;		//below the band
		bad_fixups = 0;
		while(not_optimal && (bad_fixups < MAX_BAD_FIXUPS)){
			in = !su .&& !sl;	//inside the band
			XX = selectifr(X,in);
			yy = selectifr(y,in);
		//make the globs
			if(any(su)){
				XX = XX | sumc(selectifr(X,su));
				yy = yy | sumc(selectifr(y,su));
				}
			if(any(sl)){
				XX = XX | sumc(selectifr(X,sl));
				yy = yy | sumc(selectifr(y,sl));
			}
		//reestimate with globbed sample
			b0=b;
			if(VERBOSE>=2){
				print("    Re-estimating with globbed sample of ",rows(XX)," observations\n");
			}
			b = rq_fnm(XX,yy,tau,&a);
			r = y - X*b;
		//check the globs
			qb=b-b0;
			discr=cs^2-qb'*Sinv*qb;
			
			if((bad_fixups==0 && discr<=0) || bad_fixups>0){
				su_bad = (r .<= 0) .&& su;
				sl_bad = (r .>= 0) .&& sl;
			
				if(any(su_bad ~ sl_bad)){
					su = su .&& !su_bad;
					sl = sl .&& !sl_bad;
					++bad_fixups;
					if(VERBOSE>=2){
						print("    Bad observations found!, trying again\n");
					}
				}
				else{ 
					not_optimal = FALSE;
				}
			}
			else{
				not_optimal=FALSE;
			}
		}	
	}
	return(b);
}

//**********************************************************************
// Function: rq_alltau
//
// Description:
// This function computes solutions for the primal problem (the parameter estimates)
// for all possible quantiles. The function obtains an initial solution from the
// primal-dual LP algorithm and then proceeds to iterate at the exterior of the
// constraint set to find the tau and primal solution that corresponds to every
// relevant vertex of the constraint set.
// This OX version is written by Daniel Morillo
//
// Syntax:
//		Let k be the number of tau breakpoints in the regression
//
//	Input:
//		y, n column vector of observations of the dependent variable
//		X, x by p matrix of n observations of p independent variables.
//		   The function DOES NOT automatically include a constant.
//		taus, address
//			in: nothing
//			out: k column vector containing tau breakpoints
//		bhats, address
//			in: nothing
//			out: k by p matrix of coefficients at each breakpoint.
//				 The coefficients at breakpoint s correspond to all taus
//				 from breakpoint s to breakpoint s+1 (excluding the right end)
//				 for s=1,...,k-1
//		aps, address
//			in: nothing
//			out:if DUAL_SAVEMODE=0, nothing
//			    if DUAL_SAVEMODE=1, k by n matrix of dual solutions at each breakpoint.
//			    if DUAL_SAVEMODE=2, k by 2p matrix. First p columns contain basis observation indices.
//						Second p columns contain the corresponding duals.
//		V, address
//			in: nothing
//			out: Objective function values at each tau.
//
//
//	Output:
//		-1 if there was an error
//		1 if successful
//
//**********************************************************************

rq_alltau(X,y,taus,bhats,aps,V){
	decl n,p,eps,t;
	decl starttau,tau;
	decl a,ap,oldap,bhat,bhatp,r;
	decl out,tmp1,tmp2;
	decl ag,dg,Xvar,Xvar1;
	decl tm,tm1,tm2,tm3,tm4,tm5;
	decl basic,rest;
	decl b,minima,flag;
	decl Xhinv,dn,ratio,minratio,in;
	decl ok;

	tm=timer();
	tm2=0;
	tm3=0;
	tm4=0;
	tm5=0;

	eps=1e-12;
	n=rows(y);
	p=columns(X);

	ok=0;
	while(!ok){

	taus[0]=new matrix[0][1];
	bhats[0]=new matrix[0][p];
	V[0]=new matrix[0][1];
	if(DUAL_SAVEMODE==1){
		aps[0]=new matrix[0][n];
	}
	if(DUAL_SAVEMODE==2){
		aps[0]=new matrix[0][2*p];
	}

	//get first (largest) tau

	starttau=1-(1/(n+1));

	bhat=rq_fnm(X,y,starttau,&a);

	//now "purify" the results so that we are exactly at vertex
	r=y-X*bhat;
	minima=sortc(fabs(r));
	basic=new matrix[0][1];
	for(t=1;t<=p;++t){
		basic=basic|vecindex(fabs(r).==minima[(t-1)][]);
	}
	rest=vecindex(fabs(r).>=0);
	rest=deleter(rest,basic);

	Xvar=X'*ones(n,1);
	Xhinv=X[basic][]^-1;
	bhatp=Xhinv*y[basic][];		//bhatp is purified primal solution
	r=y-X*bhatp;			//"exact" residuals

	//Now compute purified dual solution
	ap=a;
	ap[rest][]=a[rest][].<0.5 .?0 .:1;
	Xvar1=X[rest][]'*ap[rest][];
	ag=Xhinv'*(Xvar-Xvar1);
	dg=-Xhinv'*Xvar;
	ap[basic][]=ag+dg*starttau;
	r[basic][]=0;

	//start iterating

	tau=starttau;

	while(tau>eps){
		oldap=ap;
	
		//now compute new tau
		tm1=timer();

		tmp1=-ag./dg;
		tmp1=tmp1.<tau .?tmp1 .: -1e20;

		tmp2=(1-ag)./dg;
		tmp2=tmp2.<tau .?tmp2 .: -1e20;

		tau=max(tmp1~tmp2);

		//now compute which obs "is out" of the solution
		out=vecindex(tmp1~tmp2.==tau);
		flag=double(out>(p-1) ? -1 : 1);
		out=out>(p-1) ?out-p :out;
		tau=double(tau-(eps+eps/fabs(1+dg[out][]))); //go up
		out=basic[out][];

		if(tau<eps){
			tau=0;
                        taus[0]=tau|taus[0];
			V[0]=sumc(r.*(r.<0 .? tau-1 .: tau))|V[0];
                        bhats[0]=bhatp'|bhats[0];
			if(DUAL_SAVEMODE==1){
		                aps[0]=ap'|aps[0];
			}
			if(DUAL_SAVEMODE==2){
				aps[0]=(basic'~(ap[basic][])')|aps[0];
			}
			ok=1;
			break;
		}

		taus[0]=tau|taus[0];
		V[0]=sumc(r.*(r.<0 .? tau-1 .: tau))|V[0];

		bhats[0]=bhatp'|bhats[0];
		if(DUAL_SAVEMODE==1){
			aps[0]=ap'|aps[0];
		}
		if(DUAL_SAVEMODE==2){
			aps[0]=(basic'~(ap[basic][])')|aps[0];
		}

		tm4=tm4+(timer()-tm1);
		

		//find new member of basic solution
		tm1=timer();			
		dn=flag*(X*Xhinv)[][vecindex(basic.==out)];
		ratio= round((dn.==0 .? 0 .: r./dn)*1e12)/1e12; //only distinguish up to the 12th decimal
		minratio=min(ratio[rest][].>0 .? ratio[rest][] .: 1e20);
		in=vecindex(ratio.==minratio);

		if(rows(in)>1){
			//more that one obs wants to come into the basis
			//perturb the offending obs to avoid problems
			X[in][]=X[in][]+rann(rows(in),1)*1e-12;
			print("\nWarning: degeneracy detected... restarting\n");
			break;
		}

		//update basis observations
		basic=deleter(basic,out)|in;
		rest=deleter(rest,in)|out;
		tm2=tm2+(timer()-tm1);

		tm1=timer();
		Xhinv=X[basic][]^-1;
		bhatp=Xhinv*y[basic][];
		r=y-X*bhatp;

		ap[rest][]=r[rest][].<eps .?0 .:1;
		Xvar1=Xvar1+X[out][]'*ap[out][]-X[in][]'*ap[in][];
		ag=Xhinv'*(Xvar-Xvar1);
		dg=-Xhinv'*Xvar;
		ap[basic][]=ag+tau*dg;
		r[basic][]=0;
			
		tm3=tm3+(timer()-tm1);

		//check for optimality
		if(sumc(ap[deleter(rest,out)][]-oldap[deleter(rest,out)][])!=0){
			print("\nWarning: no optimality\n");
			print(sumc(ap[deleter(rest,out)][]-oldap[deleter(rest,out)][]));
			return(-1);
		}
		
	}
	}

	//tm=timer()-tm;
	//print("total time:             ",tm,"\n");
	//print(" % computing in obs:    ",tm2/tm,"\n");
	//print(" % updating variables:  ",tm3/tm,"\n");
	//print(" % computing out obs:   ",tm4/tm,"\n");
	return(0);
}

//*********************************************************************
// Function: pprog_rhs
//
// Description:
// Performs a right-hand side parametric programming step. That is, it
// computes the next incoming and outgoing observations and the boundary
// for theta. The function is desgined to be used iteratively in the computation
// of confidence intervals for quantile regression.
// This OX version is written by Daniel Morillo
//
// Syntax:
//
// 	Inputs:
//		y, n column vector of observations of the dependent variable
//		X, x by p matrix of n observations of p independent variables. No constant is included//
//		tau, scalar in (0,1) corresponding to the quantile being estimated
//		init, scalar. It tells the algorihtm if it is being called to do the first
//			iteration of the problem (1=true, 0=false)
//		bhatp, address:
//			in:	solution to the reduce problem before the step
//			out: solution to the reduced problem after the step
//		bj, address:
//			in: it contains the coefficient that we are moving from the rhs to the lhs
//			out:
//		r, address:
//			in:
//			out: latest residual vector
//		ap, address:
//			in:	dual vector before the step
//			out: dual vector after the step
//		select, address:
//			in:	positions of the variables that remain in the rhs before the step
//			out: positions of the variables that remain in the rhs after the step
//		var, scalar. It indicates the variable that we are manipulating
//		basic, address:
//			in: positions of the basis obs before the step
//			out: positions of the basis obs	after the step
//		rest, address:
//			in: positions of the non-basis obs before the step
//			out: positions of the non-basis obs after the step
//		Xhinv, address:
//			in: inverse matrix of the basis obs before the step
//			out: inverse matrix of the basis obs after the step
//		Xvar, matrix. The only part of the feasibility eqn that does not
//			change as we iterate
//		ntheta, address:
//			it conatins the maximum change in the lhs coefficients
//			that is allowed before a change of basis
//		up, scalar. It tells the algorithm if we are moving the variable up or down
//			(1=up, 0=down)
//
//	Output:
//		none
//
//*********************************************************************

pprog_rhs(const X,const y,const tau,const init, const bhatp, const bj,const r, const ap,const select, const var, const basic, const rest,const Xhinv, const Xvar,const ntheta,up){
	decl g,s,tmp,tmp1,tmp2,tmp3;
	decl in,out;
	decl newa,n,p;
	decl dir,v;
	decl d3,eps;

	n=rows(y);
	if(up==1){
		eps=1e-6;
	}
	else{
		eps=-1e-6;
	}

	if(init==1){ //first pivot: need to reduce original basis!
		p=columns(X);
		//Assign initial coefficient value

		bj[0]=bhatp[0][var][]+eps;
		bhatp[0]=deleter(bhatp[0],bj[0]);

		//obtain column of slope that corresponds to variable of interest
		dir=(Xhinv[0]')[][var];

		//now find which original basis obs would hit the boundary first
		if(up==1){
			tmp=(dir.<0 .? (-ones(p,1)+ap[0][basic[0]][])./dir .:1e10)|(dir.>0 .? ap[0][basic[0]][]./dir .:1e10);
		}
		else{
			tmp=(dir.>0 .? (ones(p,1)-ap[0][basic[0]][])./dir .:1e10)|(dir.<0 .? -ap[0][basic[0]][]./dir .:1e10);
		}
		tmp1=min(tmp);
		out=vecindex(tmp.==tmp1);
		d3=out>(p-1) ? 0 :1;
		out=out>(p-1) ? out-p :out;
		out=basic[0][out][];

		//now rearrange the indices
		basic[0]=deleter(basic[0],out);
		rest[0]=rest[0]|out;
		Xhinv[0]=X[basic[0]][select]^-1;
	
		//Now update the dual
		ap[0][out][]=d3;
		ap[0][basic[0]][]=Xhinv[0]'*((1-tau)*Xvar[select][]-X[rest[0]][select]'*ap[0][rest[0]][]);
	}

	if(init==0){
		//now note that p is reduced!
		p=rows(basic[0]);

		//first find observation that is going in
	
		s=X[rest[0]][select]*Xhinv[0]*X[basic[0]][var]-X[rest[0]][var];
		if(up==1){
			g=s.!=0 .? -r[0][rest[0]][]./s .:1e10;
			tmp=min(g.>0 .?g .:1e10);
		}
		else{
			g=s.!=0 .? -r[0][rest[0]][]./s .:-1e10;
			tmp=max(g.<0 .?g .:-1e10);
		}
		in=rest[0][vecindex(g.==tmp)][];

		ntheta[0]=tmp+eps;

		//now compute observation that is going out
		//Recall dual of existing basis obs are always between zero and one
		v=Xhinv[0]'*X[in][select]';
		if(ap[0][in][]==1){	//if dual of incoming is one then change on it must be negative!
			tmp1=(v.<0 .? -ap[0][basic[0]][]./v .:1e10)|(v.>0 .? (ones(p,1)-ap[0][basic[0]][])./v .:1e10);
		}
		else{
			//this corresponds to the incoming having dual zero (meaning change must be positive)!
			tmp1=(v.>0 .? ap[0][basic[0]][]./v .:1e10)|(v.<0 .? -(ones(p,1)-ap[0][basic[0]][])./v .:1e10);
		}
	
		newa=min(tmp1);
		if(newa>0 && newa<1){
			out=vecindex(tmp1.==newa);
			d3=out>(p-1) ? 1 :0;
			out=out>(p-1)? out-p :out;
			out=basic[0][out][];
	
			//now update all variables
		
			basic[0]=deleter(basic[0],out)|in;
			rest[0]=deleter(rest[0],in)|out;
			Xhinv[0]=X[basic[0]][select]^-1;
	
			ap[0][out][]=d3;
		}
		else{
			ap[0][in][]=1-ap[0][in][];
		}

		ap[0][basic[0]][]=Xhinv[0]'*((1-tau)*Xvar[select][]-X[rest[0]][select]'*ap[0][rest[0]][]);

		bhatp[0]=Xhinv[0]*(y[basic[0]][]-double(ntheta[0]+bj[0])*X[basic[0]][var]);
		r[0]=(y-double(ntheta[0]+bj[0])*X[][var])-X[][select]*bhatp[0];

		r[0][basic[0]][]=0;
	
	}

}

//***********************************************************************
// Function: rq_confidence
//
// Description:
// This function first computes a quantile regression using RQ and then
// calls pprog_rhs in order to compute confidence intervals for all
// coefficients. This OX version is written by Daniel Morillo
//
// Syntax:
// 
// 	Input:
//		y, n column vector of observations of the dependent variable
//		X, n by p matrix of n observations of p independent variables. No constant is included
//		tau, scalar in (0,1) determining the quantile to estimate
//		alpha, scalar in (0,1) determining the size of the confidence interval
//		conf, address
//			in: nothing
//			out: confidence interval matrix. The matrix is organized as follows:
//				 the center column contains the parameter estimates.
//				 On the left of the center column, the first column and the third columns
//				 determine the values of the parameter after and before "jumping" beyond
//				 the critical value, respectively. The second column contains an interpolated
//				 version of the left endpoint of the interval.
//				 The columns to the right of the are similarly defined.
//		glob, scalar. 1 if preprocessing should be used to obtain original solution.
//
//	Output:
//		nothing
//***********************************************************************

rq_confidence(X,y,tau,alpha,conf,glob){
	decl qf;
	decl Xhinv,Xvar,bhat,r;
	decl basic,rest,ap,bhatp;
	decl basicold,restold,apold,bhatpold,Xhinvold,rold;
	decl ag,dg;
	decl minima,t,k,p,n;
	decl ntheta;
	decl stop,init,up;
	decl select,var;
	decl bj,Tn,Bn,Sn,qns,As;
	decl sl,pa,pb;

	n=rows(X);
	p=columns(X);
	
	conf[0]=new matrix[7][p];
	
	if(VERBOSE>=1){
		print("\nStarting Rank confidence intervals procedure\n");
	}

	//first obtain solution to original problem.
	if(glob==1){
		bhat=RQ(X,y,tau);
	}
	else{
		bhat=rq_fnm(X,y,tau,&ap);
	}
	
	//now "purify" the results so that we are exactly at vertex
	r=y-X*bhat;
	minima=sortc(fabs(r));
	basic=new matrix[0][1];
	for(t=1;t<=p;++t){
		basic=basic|vecindex(fabs(r).==minima[(t-1)][]);
	}
	rest=vecindex(fabs(r).>=0);
	rest=deleter(rest,basic);

	Xvar=X'*ones(n,1);
	Xhinv=X[basic][]^-1;
	bhatp=Xhinv*y[basic][];		//bhatp is purified primal solution
	r=y-X*bhatp;			//"exact" residuals
	r[basic][]=0;			//make sure zeros are really zeros
	
	//Now compute purified dual solution
	ap=new matrix[n][1];
	ap[rest][]=r[rest][].<0 .?0 .:1;
	ap[basic][]=Xhinv'*((1-tau)*Xvar-X[rest][]'*ap[rest][]);

	
	//Save original values
	Xhinvold=Xhinv;
	apold=ap;
	bhatpold=bhatp;
	basicold=basic;
	restold=rest;
	rold=r;

	//Now start iterating
	for(var=0;var<p;++var){
		//print operating info
		if(VERBOSE>=2){
				print("\n  Now doing PP on variable ",var+1,"\n");
		}
		
		//Now define variables not to be tested
		select=range(0,(p-1))';
		select=deleter(select,var);
	
		//Now start computing test statistic
		ols2c(X[][var],X[][select],&qf);
		//qns=(1/n)*X[][var]'*(diag(ones(n,1))*X[][var]-X[][select]*qf);
		qns=(1/n)*X[][var]'*(X[][var]-X[][select]*qf);

		As=tau*(1-tau);
		
		//now iterate for up and down directions
		for(up=0;up<=1;++up){
			if(VERBOSE>=2){
				if(up==0){
					print("    Going up\n");
				}
				else{
					print("    Going down\n");
				}
			}
			//Restore original values
			Xhinv=Xhinvold;
			ap=apold;
			bhatp=bhatpold;
			basic=basicold;
			rest=restold;
			r=rold;

			//Now loop
			ntheta=0;
			stop=0;
			init=1;
			while(stop==0){
				pprog_rhs(X,y,tau,init,&bhatp,&bj,&r,&ap,select,var,&basic,&rest,&Xhinv, Xvar,&ntheta,up);
				bj=bj+ntheta;
				Bn=ap-(1-tau);
				Sn=(1/sqrt(n))*X[][var]'*Bn;
				Tn=Sn/sqrt(As*qns);
				if(probn(fabs(Tn))>(1-alpha/2)){
					pa=probn(fabs(Tn));
					stop=1;
					if(up==1){
						conf[0][6][var]=bj;
						//interpolate
						sl=(pa-pb)/(conf[0][6][var]-conf[0][4][var]);
						conf[0][5][var]=((1-alpha/2)-pb)/sl+conf[0][4][var];
					}
					else{
						conf[0][0][var]=bj;
						//interpolate
						sl=(pa-pb)/(conf[0][0][var]-conf[0][2][var]);
						conf[0][1][var]=((1-alpha/2)-pb)/sl+conf[0][2][var];
					}

				}
				init=0;
				if(stop==0){
					pb=probn(fabs(Tn));
					if(up==1){
						conf[0][4][var]=bj;
					}
					else{
						conf[0][2][var]=bj;
					}
				}
			}
		}
	}
	conf[0][3][]=bhatpold';
	conf[0]=conf[0]';
	return(0);
}

//******************************************************************************
// Function: kerncv
//
// Description:
// This function implements the kernel method for computing variance-covariance
// matrices in quantile regression, as described in Buchinsky (1995) and originally
// proposed by Powell (1986) in the context of censored regression quantiles.
// This OX version is written by Daniel Morillo
//
// Syntax:
//
//	Input:
//		y, n column vector of observations of the dependent variable
//		X, x by p matrix of n observations of p independent variables. No constant is included
//		bhat, a p column vector of parameter estimates
//		b, positive scalar determining the exponent of the bandwidth
//
//	Output:
//		A p by p variance covariance matrix
//
//	Options:
//
//		BW_EXP: This variable controls the exponent of the bandwidth.
//
// 		BW_CONS: This variable controls the constant in the bandwidth. If set to 0, the function will
//			automatically choose this constant to be equal to the scale of the reponse variable.
//			If not zero, the constant will be set equal to the value of this variable.
//
//		KERNEL: At this point only a norma kernel has been implemented.
//
//******************************************************************************

kerncv(X,y,tau,bhat){
	decl V,n,r,h,c;

	r=y-X*bhat;
	n=rows(X);
	if(BW_CONS==0){
		c=double(quantilec(r,0.75)-quantilec(r,0.25));
	}
	else{
		c=BW_CONS;
	}

	h=(n^-BW_EXP)*c;

	if(KERNEL=="normal"){
		V=((1/(n*h))*(densn(r/h).*X)'*X)^-1;
	}
	
	V=(1/n)*tau*(1-tau)*(V*((1/n)*X'*X)*V);

	return(V);
}

//******************************************************************************
// Function: interp
//
// Description:
// This function will compute interpolated values of the variable y (which have been
// observed at values x) at values otherx. The interpolation is bounded above at the
// values of othery.
//
// Syntax:
//
//
// 	Input:
//		x, n1 column vector, sorted in ascending order
//		y, n1 column vector, sorted according to x
//		otherx, n2 column vector, sorted in ascending order
//		othery, n2 column vector, sorted according to otherx
//		
//	Output:
//		A n2 column vector corresponding to the interpolated values.
//******************************************************************************
interp(const x, const y, const otherx, const othery){
	decl newy;
	decl i,n,tmp;
	decl yb,ya,xb,xa;

	n=rows(otherx);
	newy=new matrix[n][1];

	for(i=0;i<n;++i){
		tmp=sortbyc(y~x~fabs(x-otherx[i]),2);
		yb=tmp[0][0];
		ya=tmp[1][0];
		xb=tmp[0][1];
		xa=tmp[1][1];
		
		//perform interpolation
               	newy[i]=yb+(otherx[i]-xb)*((ya-yb)/(xa-xb));
				if(newy[i]>othery[i]){
					newy[i]=othery[i];
				}
	}

	return(newy);
}
	

//**********************************************************************************************************
// function: rq_inference
//
// Description:
// This function implements some of the inference procedures for quantile regression
// found in Koenker and Machado (1999). The code is based on the "R5" S function
// written by Roger Koenker. This OX version is written by Daniel Morillo
//
// Syntax:
//
//	Input:
//		X, n by p matrix of n observations of p variables. A constant must be included
//		y, n column vector of dependent variable observations
//		tau, scalar in (0,1) corresponding to the quantile of interest. if set to -1 the entire
//			tau process will be used.
//		useXid, p1 column vector corresponding to the indices of the p1 variables in the restricted model.
//				Since X includes a constant, setting this to zero picks up only the constant
//		previous, has two formats:
//			if tau=-1 then this is a matrix containing in the first column the tau breakpoints and
//				in the second and third columns the primal solutions and thevalues of the objective
//				function at each breakpoint, respectively. Thus, in order to use this function, one
//				must first estimate the unrestricted model.
//			if tau!=-1 then this simply contains the parameter estimates (in column vector form) at tau
//		rtaus, address:
//			in: nothing
//			out: if tau=-1 a vector containing the taus breakpoints of the restricted model.
//		R1, address
//			in: 1 if R1 is desired, 0 if not;
//			out: if tau!=-1 the value of R1 at tau
//			     if tau=-1 a two column matrix with the values of the tau breakpoints for the restricted model
//				in the first column and the corresponding values of R1 in the second column.
//		Lk
//			in: 1 if the Likelihood ratio statistics are desired, o if not;
//			out: if tau!=-1 a 2 dimensional row vector with the Likelihood statistics
//				 if tau=-1 a 2 column matrix with each column containing the Likelihood statistic at each tau breakpoint
//					of the restricted model.
//	Output:
//		nothing
//**********************************************************************************************************
rq_inference(const X,const y,const tau,const useXid,const previous,const rtaus,const R1, const Lk){
	decl Vhat,Vhat1,Vtilde;
	decl utaus,ubhats;
	decl rbhats,raps;
	decl n,n1,p,i;
	decl r;
	decl L,s;

	n=rows(y);
	p=columns(X);
	
	if(tau==-1){
		utaus=previous[][0]; //retrieve tau breakpoints for unrestricted model
		ubhats=previous[][1:p]; //retrieve primal solutions for unrestricted model
		Vhat=previous[][p+1];  //retrieve Obj. fn. value for unsrestricted model
	
		//Now evaluate the restricted model
		rq_alltau(X[][useXid],y,rtaus,&rbhats,&raps,&Vtilde);

		n1=rows(Vtilde);
		
		//now interpolate Vhat to obtain Vhat at rtaus

		Vhat1=interp(utaus,Vhat,rtaus[0],Vtilde);

		if(R1[0]==1){ //Compute R1
			R1[0]=1-Vhat1./Vtilde; //NOTE; first value will be NaN since objective at tau=0 is zero for both constrained and unconstarined models
			R1[0]=R1[0][1:rows(rtaus[0])-1][];
		}
		if(Lk[0]==1){ //Compute LR test using Laplace densities
			Lk[0]=new matrix[n1-1][2];
			for(i=1;i<n1;++i){
				//compute the sparsity
				s=sparsity(y-X*(ubhats[i][])',rtaus[0][i][]);
				//Compute the LR statistics
				Lk[0][i-1][0]=2*(Vtilde[i][]-Vhat1[i][])/((rtaus[0][i][]*(1-rtaus[0][i][]))*s);
				Lk[0][i-1][1]=(2*Vhat1[i][]/(rtaus[0][i][]*(1-rtaus[0][i][])*s))*log(Vtilde[i][]/Vhat1[i][]);
			}
		}
	}
	else{
		r=y-X*previous;
		Vhat=sumc(r.*(r.<0 .? tau-1 .: tau));

		//Now evaluate the restricted model
        rbhats=rq_fnm(X[][useXid],y,tau,&raps);
		r=y-rbhats*(X[][useXid]);
	
		Vtilde=sumc(r.*(r.<0 .? tau-1 .: tau));

		if(R1[0]==1){
			R1[0]=1-Vhat/Vtilde;
		}
		if(Lk[0]==1){ //Compute LR test using Laplace densities
			//compute the sparsity
			s=sparsity(y-X*previous,tau);
			Lk[0][0]=2*(Vtilde-Vhat)/(tau*(1-tau)*s);
			Lk[0][1]=(2*Vhat/(tau*(1-tau)*s))*log(Vtilde/Vhat);
		}
	}
	
}


//******************************************************************************
// Function: rq
//
// Description:
// This is just a wrapper function to allow the user to choose among several different
// estimation options
//
// Syntax:
//
//	Input:
//		y, n column vector of observations of the dependent variable
//		X, x by p matrix of n observations of p independent variables. A constant is automatically included
//		tau, scalar in (0,1) corresponding ot the quantile to estimate. If the entire tau
//		     process is required, set this to -1
//
//	Output:
//		A 7 element array if successful. -1 if there's an error.
//		The indices of the array denote:
//		0: if RCI=0, a p column vector of parameter estimates.
//		   if RCI=1, a p by 3 matrix containg the estimates in the center column and the lower and
//		             upper bounds of the ci's in the first and third columns respectively.
//		   if tau=-1, nothing
//
//		1: if tau!=-1, nothing
//		   if tau=-1, a column vector containing all k breakpoints of the tau proceess
//
//		2: if tau!=-1, nothing
//		   if tau=-1, a k by p matrix containing in each row the parameter estimates
//				corresponding to each tau in the tau process.
//
//		3: if tau!=-1, nothing
//		   if tau=-1, a k by n matrix containing in each row the dual solution for the
//				corresponding tau in the tau process. See rq_alltau for the format
//				in which the dual solution is reported.
//
//		4: if tau!=-1, nothing
//		   if tau=-1, a k by n matrix containing the objective function value at each tau breakpoint
//		
//		5: if KCV=1, a p by p variance covariance matrix for the paramater estimates
//		   if KCV=0 or tau=-1, nothing
//
//		6: if REPORT_R1=0, nothing
//		   if REPORT_R1=1, there are two formats:
//			if tau=-1, a two-column matrix with the values of the tau breakpoints of the restricted model
//					in the first column and the corresponding R1 values in the second column.
//			if tau!=-1, a scalar corresponding to the value of R1 at tau.
//
//	Options:
//		GLOB, scalar. 1 for preprocessing, 0 for no preprocessing. Preprocessing is only
//			  an option if tau!=-1
//		RCI, scalar. 1 for rank confidence itervals, 0 for no rank confidence intervals
//		ALPHA, scalar in (0,1) determining the size of the confidence interval
//		KCV, scalar. 1 for kernel var-covar matrix, 0 for no var-covar matrix
//********************************************************************************
rq(X,y,tau){
	decl res,result,tmp1,tmp2,tmp3,tmp4,n,bhat;

	n=rows(X);

	result=new array[7];

	if(VERBOSE>=1){
		print("\nStarting estimation procedure\n");
		print("Options:\n");
		if(tau<0){
			print("  Computing entire tau proces.Dual save mode is ",DUAL_SAVEMODE,"\n");
			if(REPORT_R1==1){
				print("  Computing entire R1 process\n");
			}
		}
		else{
			print("  Computing single tau=",tau,"\n");
			if(GLOB==1){
				print("  Globbing procedures will be used\n");
			}
			if(RCI==1){
				print("  Rank confidence Intervals will be computed\n");
			}
			if(KCV==1){
				print("  Kernel Var-Covar matrix will be computed\n");
			}
			if(REPORT_R1){
				print("  R1 will be computed\n");
			}
		}
		print("  Verbose mode is ",VERBOSE,"\n");
	}
	
	if(tau==-1){
		res=rq_alltau(ones(n,1)~X,y,&tmp1,&tmp2,&tmp3,&tmp4);
		if(res==-1){
			result=-1;
		}
		else{
			result[1]=tmp1; //tau breakpoints
			result[2]=tmp2; //Primal solutions
			result[3]=tmp3;	//Dual solutions
			result[4]=tmp4; //Objective function values
		}
		tmp1=result[1]~result[2]~result[4]; //for use in calling rq_inference below
	}
	else{
		if(RCI==1){
			rq_confidence(ones(n,1)~X,y,tau,ALPHA,&tmp1,GLOB);
			result[0]=tmp1[][1]~tmp1[][3]~tmp1[][5];
			bhat=tmp1[][3];
		}
		else{
			if(GLOB==1){
				result[0]=RQ(ones(n,1)~X,y,tau);
				bhat=result[0];
			}
			else{
				result[0]=rq_fnm(ones(n,1)~X,y,tau,&tmp1);
				bhat=result[0];
			}
		}
		if(KCV==1){
			result[5]=kerncv(ones(n,1)~X,y,tau,bhat);
		}

		tmp1=bhat; // for use in calling rq_inference below
	}

	if(REPORT_R1==1){
	   	tmp3=1;
		rq_inference(ones(n,1)~X,y,tau,0,tmp1,&tmp2,&tmp3,&tmp4);
		if(tau==-1){
			result[6]=tmp2[1:rows(tmp2)-1]~tmp3;
		}
		else{
			result[6]=tmp3;
		}
	}

	return(result);
}
