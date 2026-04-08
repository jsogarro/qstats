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
