/ distributions.q — Probability Distributions for qstats
/ qstats: Statistical Testing & Diagnostics for Q/kdb+
/ License: MIT
/ Namespace: .dist

/ @desc Normal probability density function
/ @param xx:float — input value(s)
/ @param mu:float — mean
/ @param sigma:float — standard deviation
/ @return float — PDF value(s)
.dist.dnorm:{[xx;mu;sigma]
  z:(xx-mu)%sigma;
  coef:1f%(sigma*sqrt 2*.special.PI);
  coef*exp neg 0.5*z*z
 };

/ @desc Normal cumulative distribution function via Abramowitz & Stegun
/ @param xx:float — input value(s)
/ @param mu:float — mean
/ @param sigma:float — standard deviation
/ @return float — CDF value(s)
.dist.pnorm:{[xx;mu;sigma]
  / @desc Normal CDF via Abramowitz & Stegun rational approximation
  / @param x - value (scalar or vector), mu - mean, sigma - std dev
  / @return cumulative probability
  zz:(xx-mu)%sigma;
  az:abs zz;
  tt:1%1+0.2316419*az;
  dd:0.3989422804*exp neg 0.5*az*az;
  pp:dd*tt*(0.3193815+tt*(-0.3565638+tt*(1.781478+tt*(-1.821256+tt*1.330274))));
  ?[zz<0;pp;1-pp]
 };

/ @desc Normal quantile function via Acklam's algorithm (scalar)
/ @param p:float — probability (0 < p < 1)
/ @param mu:float — mean
/ @param sigma:float — standard deviation
/ @return float — quantile value
.dist.qnorm_scalar:{[p;mu;sigma]
  / Acklam coefficients
  a1:-3.969683028665376e+01;
  a2:2.209460984245205e+02;
  a3:-2.759285104469687e+02;
  a4:1.383577518672690e+02;
  a5:-3.066479806614716e+01;
  a6:2.506628277459239e+00;
  /
  b1:-5.447609879822406e+01;
  b2:1.615858368580409e+02;
  b3:-1.556989798598866e+02;
  b4:6.680131188771972e+01;
  b5:-1.328068155288572e+01;
  /
  c1:-7.784894002430293e-03;
  c2:-3.223964580411365e-01;
  c3:-2.400758277161838e+00;
  c4:-2.549732539343734e+00;
  c5:4.374664141464968e+00;
  c6:2.938163982698783e+00;
  /
  d1:7.784695709041462e-03;
  d2:3.224671290700398e-01;
  d3:2.445134137142996e+00;
  d4:3.754408661907416e+00;
  /
  plo:0.02425;
  phi:1-plo;
  /
  / For scalar p: use conditional logic
  / Central region
  $[(p>=plo) & p<=phi;
    [qq:p-0.5; rr:qq*qq;
     numer:(a1*rr)+a2; numer:(numer*rr)+a3; numer:(numer*rr)+a4; numer:(numer*rr)+a5; numer:(numer*rr)+a6; numer:numer*qq;
     denom:(b1*rr)+b2; denom:(denom*rr)+b3; denom:(denom*rr)+b4; denom:(denom*rr)+b5; denom:(denom*rr)+1;
     mu+(sigma*numer)%denom];
    p<plo;
    [rr:sqrt neg 2*log p;
     numer:(c1*rr)+c2; numer:(numer*rr)+c3; numer:(numer*rr)+c4; numer:(numer*rr)+c5; numer:(numer*rr)+c6;
     denom:(d1*rr)+d2; denom:(denom*rr)+d3; denom:(denom*rr)+d4; denom:(denom*rr)+1;
     mu+(sigma*numer)%denom];
    / p > phi (upper tail)
    [rr:sqrt neg 2*log[1-p];
     numer:(c1*rr)+c2; numer:(numer*rr)+c3; numer:(numer*rr)+c4; numer:(numer*rr)+c5; numer:(numer*rr)+c6;
     denom:(d1*rr)+d2; denom:(denom*rr)+d3; denom:(denom*rr)+d4; denom:(denom*rr)+1;
     mu+(sigma*neg numer)%denom]
  ]
 };
/ @desc Normal quantile function via Acklam's algorithm
/ @param pp:float — probability (0 < p < 1) - scalar or vector
/ @param mu:float — mean
/ @param sigma:float — standard deviation
/ @return float — quantile value(s)
.dist.qnorm:{[pp;mu;sigma]
  .dist.qnorm_scalar[;mu;sigma] each pp
 };

/ @desc Normal random variates via Box-Muller transform
/ @param nn:long — number of samples
/ @param mu:float — mean
/ @param sigma:float — standard deviation
/ @return float[] — random samples
.dist.rnorm:{[nn;mu;sigma]
  u1:nn?1f;
  u2:nn?1f;
  r:sqrt neg 2f*log u1;
  theta:2f*.special.PI*u2;
  zz:r*cos theta;
  mu+sigma*zz
 };

/=============================================================================
/ CHI-SQUARED DISTRIBUTION
/=============================================================================

/ @desc Chi-squared probability density function
/ @param xx:float — input value(s) (>= 0)
/ @param df:float — degrees of freedom (> 0)
/ @return float — PDF value(s)
.dist.dchisq:{[xx;df]
  a:df%2f;
  / Compute log PDF for positive x
  t1:(a-1f)*log xx;
  t2:xx%2f;
  t3:a*log 2f;
  t4:.special.lgamma a;
  lnf:((t1-t2)-t3)-t4;
  / Handle x=0 edge case
  edge_val:$[df>2f;0f;$[df=2f;0.5;0w]];
  ?[xx<=0f;edge_val;exp lnf]
 };

/ @desc Chi-squared cumulative distribution function
/ @param xx:float — input value(s) (>= 0)
/ @param df:float — degrees of freedom (> 0)
/ @return float — CDF value(s)
.dist.pchisq:{[xx;df]
  a:df%2f;
  .special.gammainc[a;] each xx%2f
 };

/ @desc Chi-squared quantile function via Newton-Raphson
/ @param pp:float — probability (0 < p < 1)
/ @param df:float — degrees of freedom (> 0)
/ @return float — quantile value(s)
.dist.qchisq_scalar:{[p;df]
  / Special case: df=1 has exact formula
  $[df<1.5;
    [z1:.dist.qnorm_scalar[(1f+p)%2f;0f;1f]; z1*z1];
    / General case: Wilson-Hilferty approximation
    [z:.dist.qnorm_scalar[p;0f;1f];
     wh_inner:(1f-(2f%(9f*df)))+(z*sqrt 2f%(9f*df));
     / Use simpler guess if Wilson-Hilferty gives negative value
     x:$[wh_inner>0f; df*xexp[wh_inner;3]; $[p<0.1; df*p; df]];
     / Newton-Raphson iterations
     i:0; maxiter:100; tol:1e-10;
     while[(i<maxiter) and tol<abs fx:.dist.pchisq[x;df]-p;
       dpdf:.dist.dchisq[x;df];
       / Guard against zero/tiny derivative
       x:$[dpdf>1e-200; x-(fx%dpdf); x*0.5];
       i+:1];
     x]
  ]
 };
.dist.qchisq:{[pp;df] .dist.qchisq_scalar[;df] each pp};

/ @desc Chi-squared random variates via sum of squared normals
/ @param nn:long — number of samples
/ @param df:float — degrees of freedom (> 0)
/ @return float[] — random samples
.dist.rchisq:{[nn;df]
  / Generate df × nn standard normals, square and sum across df
  z:.dist.rnorm[df*nn;0f;1f];
  z2:z*z;
  sum each (df;nn)#z2
 };

/=============================================================================
/ STUDENT'S T-DISTRIBUTION
/=============================================================================

/ @desc Student's t probability density function
/ @param xx:float — input value(s)
/ @param df:float — degrees of freedom (> 0)
/ @return float — PDF value(s)
.dist.dt:{[xx;df]
  / Use log-gamma to avoid overflow
  lg1:exp .special.lgamma (df+1f)%2f;
  lg2:exp .special.lgamma df%2f;
  numer:lg1;
  denom:(sqrt df*acos -1f)*lg2;
  coef:numer%denom;
  coef*xexp[1f+(xx*xx)%df;(neg df+1f)%2f]
 };

/ @desc Student's t cumulative distribution function
/ @param xx:float — input value(s)
/ @param df:float — degrees of freedom (> 0)
/ @return float — CDF value(s)
.dist.pt:{[xx;df]
  z:df%(df+xx*xx);
  a:df%2f;
  b:0.5;
  betacdf:0.5*.special.betainc[;a;b] each z;
  ?[xx<0;betacdf;1f-betacdf]
 };

/ @desc Student's t quantile function via Newton-Raphson
/ @param pp:float — probability (0 < p < 1)
/ @param df:float — degrees of freedom (> 0)
/ @return float — quantile value(s)
.dist.qt_scalar:{[p;df]
  / Initial guess: use normal quantile
  x:.dist.qnorm_scalar[p;0f;1f];
  / Newton-Raphson iterations
  i:0; maxiter:50; tol:1e-10;
  while[(i<maxiter) and tol<abs fx:.dist.pt[x;df]-p;
    x:x-(fx%.dist.dt[x;df]);
    i+:1];
  x
 };
.dist.qt:{[pp;df] .dist.qt_scalar[;df] each pp};

/ @desc Student's t random variates via ratio method
/ @param nn:long — number of samples
/ @param df:float — degrees of freedom (> 0)
/ @return float[] — random samples
.dist.rt:{[nn;df]
  z:.dist.rnorm[nn;0f;1f];       / Standard normal
  v:.dist.rchisq[nn;df];         / Chi-squared
  z%sqrt v%df
 };

/=============================================================================
/ F-DISTRIBUTION
/=============================================================================

/ @desc F-distribution probability density function
/ @param xx:float — input value(s) (>= 0)
/ @param df1:float — numerator degrees of freedom (> 0)
/ @param df2:float — denominator degrees of freedom (> 0)
/ @return float — PDF value(s)
.dist.df:{[xx;df1;df2]
  a:df1%2f;
  b:df2%2f;
  lnbeta:(.special.lgamma a)+(.special.lgamma b)-.special.lgamma a+b;
  / F PDF log form: a*log(d1) + b*log(d2) + (a-1)*log(x) - (a+b)*log(d1*x+d2) - log(Beta(a,b))
  t1:a*log df1;
  t2:b*log df2;
  t3:(a-1f)*log xx;
  t4:(a+b)*log (df1*xx)+df2;
  t5:lnbeta;
  lnf:((t1+t2)+t3)-(t4+t5);
  exp lnf
 };

/ @desc F-distribution cumulative distribution function
/ @param xx:float — input value(s) (>= 0)
/ @param df1:float — numerator degrees of freedom (> 0)
/ @param df2:float — denominator degrees of freedom (> 0)
/ @return float — CDF value(s)
.dist.pf:{[xx;df1;df2]
  z:(df1*xx)%(df1*xx)+df2;
  a:df1%2f;
  b:df2%2f;
  .special.betainc[;a;b] each z
 };

/ @desc F-distribution quantile function via Newton-Raphson
/ @param pp:float — probability (0 < p < 1)
/ @param df1:float — numerator degrees of freedom (> 0)
/ @param df2:float — denominator degrees of freedom (> 0)
/ @return float — quantile value(s)
.dist.qf_scalar:{[p;df1;df2]
  / Initial guess: ratio of chi-squared quantiles
  x:(.dist.qchisq_scalar[p;df1]%df1)%.dist.qchisq_scalar[1f-p;df2]%df2;
  / Newton-Raphson iterations
  i:0; maxiter:50; tol:1e-10;
  while[(i<maxiter) and tol<abs fx:.dist.pf[x;df1;df2]-p;
    x:x-(fx%.dist.df[x;df1;df2]);
    i+:1];
  x
 };
.dist.qf:{[pp;df1;df2] .dist.qf_scalar[;df1;df2] each pp};

/ @desc F-distribution random variates via ratio of chi-squareds
/ @param nn:long — number of samples
/ @param df1:float — numerator degrees of freedom (> 0)
/ @param df2:float — denominator degrees of freedom (> 0)
/ @return float[] — random samples
.dist.rf:{[nn;df1;df2]
  v1:.dist.rchisq[nn;df1];
  v2:.dist.rchisq[nn;df2];
  (v1%df1)%(v2%df2)
 };

/=============================================================================
/ UNIFORM DISTRIBUTION
/=============================================================================

/ @desc Uniform probability density function
/ @param xx:float — input value(s)
/ @param a:float — lower bound (default 0)
/ @param b:float — upper bound (default 1)
/ @return float — PDF value(s)
.dist.dunif:{[xx;a;b]
  ?[(xx>=a)&xx<=b;1f%(b-a);0f]
 };

/ @desc Uniform cumulative distribution function
/ @param xx:float — input value(s)
/ @param a:float — lower bound (default 0)
/ @param b:float — upper bound (default 1)
/ @return float — CDF value(s)
.dist.punif:{[xx;a;b]
  ?[xx<a;0f;?[xx>b;1f;(xx-a)%(b-a)]]
 };

/ @desc Uniform quantile function
/ @param pp:float — probability (0 <= p <= 1)
/ @param a:float — lower bound (default 0)
/ @param b:float — upper bound (default 1)
/ @return float — quantile value(s)
.dist.qunif:{[pp;a;b]
  a+pp*(b-a)
 };

/ @desc Uniform random variates
/ @param nn:long — number of samples
/ @param a:float — lower bound (default 0)
/ @param b:float — upper bound (default 1)
/ @return float[] — random samples
.dist.runif:{[nn;a;b]
  a+(b-a)*nn?1f
 };
