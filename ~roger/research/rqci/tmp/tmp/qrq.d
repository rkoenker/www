.BG
.FN qrq
.TL
Linearized Quantile Estimation
.DN
Compute linearized quantiles from rq data structure.
.CS
qrq(s, a)
.RA
.AG s
data structure returned by the quantile regression function rq with t<0 or t>1.
.AG a
the vector of quantiles for which the corresponding linearized quantiles are to be computed.
.RT
a vector of the linearized quantiles corresponding to vector a,
as interpolated from the second row of s$sol.
.SA SEE ALSO
rq and  trq  for further detail.
.EX
z_qrq(rq(x,y),a)       #assigns z the linearized quantiles 
                       #corresponding to vector a.
.KW 
.WR
