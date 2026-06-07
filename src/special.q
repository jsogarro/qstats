/ special.q — Special Functions for qstats
/ qstats: Statistical Testing & Diagnostics for Q/kdb+
/ License: MIT
/ Namespace: .special

.special.PI:acos -1f;

/ @desc Log-gamma function via Lanczos approximation (g=7, n=9)
/ @param z:float — input value(s)
/ @return float — ln(Gamma(z))
.special.lgamma:{[z_input]
  / Lanczos coefficients for g=7, n=9
  c0:0.99999999999980993;
  c1:676.5203681218851;
  c2:-1259.1392167224028;
  c3:771.32342877765313;
  c4:-176.61502916214059;
  c5:12.507343278686905;
  c6:-0.13857109526572012;
  c7:9.9843695780195716e-6;
  c8:1.5056327351493116e-7;
  g:7f;
  /
  / Handle reflection for z < 0.5
  / For very small positive z (<0.1), use recurrence first: lgamma(z) = lgamma(z+1) - ln(z)
  use_recurrence:(z_input>0)&z_input<0.1;
  zz_adj:?[use_recurrence;z_input+1;z_input];
  recur_corr:?[use_recurrence;neg log z_input;0f];
  /
  use_reflection:zz_adj<0.5;
  zz:?[use_reflection;1f-zz_adj;zz_adj];
  /
  / Compute Ag(z) = c0 + sum(ck/(z+k-1))
  ag:c0;
  ag+:c1%(zz+0f);
  ag+:c2%(zz+1f);
  ag+:c3%(zz+2f);
  ag+:c4%(zz+3f);
  ag+:c5%(zz+4f);
  ag+:c6%(zz+5f);
  ag+:c7%(zz+6f);
  ag+:c8%(zz+7f);
  /
  / Compute lgamma via Lanczos formula
  term1:0.5*log 2*.special.PI;
  term2:(zz-0.5)*log zz+g-0.5;
  term3:neg zz+g-0.5;
  term4:log ag;
  result:term1+term2+term3+term4;
  /
  / Apply reflection formula for z < 0.5
  / lgamma(z) = ln(pi) - ln(|sin(pi*z)|) - lgamma(1-z)
  / Note: q's `-` is right-associative, so a - b - c parses as a - (b - c).
  / Parenthesise the first subtraction to get the intended (a - b) - c.
  reflection_term:((log .special.PI) - (log abs sin .special.PI*zz_adj)) - result;
  final:?[use_reflection;reflection_term;result];
  / Apply recurrence correction
  final+recur_corr
 };

/ @desc Regularized incomplete beta function I_x(a,b) via Lentz continued fraction.
/ Helper for .special.betainc. Use the dispatcher (.special.betainc) directly
/ unless you have a specific reason; this helper assumes the caller has already
/ applied symmetry so the relevant x is in the smaller half.
/ @param xx:float — upper limit (0 < x < 1)
/ @param aa:float — shape parameter a > 0
/ @param bb:float — shape parameter b > 0
/ @return float — I_x(a,b)
.special.betainc_cf:{[xx;aa;bb]
  / Front factor: x^a * (1-x)^b / B(a,b)
  lnfront:(aa*log[xx])+(bb*log[1-xx]);
  lnfront:lnfront+.special.lgamma[aa+bb];
  lnfront:lnfront-.special.lgamma[aa];
  lnfront:lnfront-.special.lgamma[bb];
  front:exp lnfront;
  / Continued fraction via Lentz's method
  / The CF is: 1/(1+ d1/(1+ d2/(1+ d3/(1+ ...))))
  / where d_{2m} = m*(b-m)*x / ((a+2m-1)*(a+2m))
  /       d_{2m+1} = -(a+m)*(a+b+m)*x / ((a+2m)*(a+2m+1))
  tiny:1e-30;
  ff:1f;
  cc:1f;
  dd:0f;
  ii:1;
  while[ii<=200;
    mm:ii div 2;
    dm:$[0=ii mod 2;
      / Even: d_{2m} = m*(b-m)*x / ((a+2m-1)*(a+2m))
      (mm*(bb-mm)*xx) % ((aa+ii-1)*(aa+ii));
      / Odd: d_{2m+1} = -(a+m)*(a+b+m)*x / ((a+2m)*(a+2m+1))
      (neg (aa+mm)*(aa+bb+mm)*xx) % ((aa+ii-1)*(aa+ii))
    ];
    dd:1+dm*dd;
    if[(abs dd)<tiny; dd:tiny];
    dd:1%dd;
    cc:1+dm%cc;
    if[(abs cc)<tiny; cc:tiny];
    delta:cc*dd;
    ff:ff*delta;
    if[1e-12>abs delta-1; ii:201];
    ii+:1;
  ];
  front%(aa*ff)
 };

/ @desc Regularized incomplete beta via power series. Numerically stable for
/ small parameters and x in the small half of [0, 1] where the Lentz CF
/ collapses (mixing positive/negative terms with poor cancellation).
/ Identity: I_x(a,b) = (x^a (1-x)^b / (a B(a,b))) * 2F1(a+b, 1; a+1; x)
/ Hypergeometric series: sum_{k=0}^inf (a+b)_k/(a+1)_k * x^k with term
/ recurrence t_k = t_{k-1} * (a+b+k-1)/(a+k) * x.
/ Geometric convergence ratio approaches x, so x <= 0.5 ensures convergence
/ within 50 iters at 1e-15 precision.
.special.betainc_series:{[xx;aa;bb]
  / Front factor: x^a * (1-x)^b / (a * B(a,b)).
  lnfront:(aa*log[xx])+(bb*log[1f-xx]);
  lnfront+:.special.lgamma[aa+bb];
  lnfront-:log[aa];
  lnfront-:.special.lgamma[aa];
  lnfront-:.special.lgamma[bb];
  front:exp lnfront;
  s:1f;
  term:1f;
  k:1;
  while[k<200;
    term*:((aa+bb+k-1f)%(aa+k))*xx;
    s+:term;
    if[(abs term)<1e-15*abs s; k:200];
    k+:1;
  ];
  front*s
 };

/ @desc Regularized incomplete beta function I_x(a,b).
/ Dispatcher: applies the (1-x, b, a) symmetry to keep x in the smaller half,
/ then picks the power series (numerically stable for small parameters and
/ small x) or the Lentz CF based on Cephes' standard gate (b*x <= 1 AND
/ x <= 0.5). The earlier impl used CF unconditionally, which lost precision
/ for a, b <= 0.5 -- producing non-monotone outputs visible at a=b=0.1.
/ @param xx:float — upper limit (0 <= x <= 1)
/ @param aa:float — shape parameter a > 0
/ @param bb:float — shape parameter b > 0
/ @return float — I_x(a,b)
.special.betainc:{[xx;aa;bb]
  if[xx<=0f; :0f];
  if[xx>=1f; :1f];
  if[xx>(aa+1)%(aa+bb+2); :1f-.special.betainc[1f-xx;bb;aa]];
  / Note: ((bb*xx)<=1f), not (bb*xx<=1f) — q right-to-left would otherwise
  / parse the latter as `bb * (xx<=1f)`, multiplying by a bool.
  $[((bb*xx)<=1f) and xx<=0.5;
    .special.betainc_series[xx;aa;bb];
    .special.betainc_cf[xx;aa;bb]]
 };

/ @desc Regularized incomplete gamma P(a,x) - Series expansion helper
/ @param aa:float — shape parameter a > 0
/ @param xx:float — upper limit of integration x >= 0
/ @return float — P(a,x) via series
.special.gammainc_series:{[aa;xx]
  / Series: P(a,x) = exp(-x + a*ln(x) - lgamma(a)) * sum(x^n / prod(a+k))
  lnpre:(neg xx)+(aa*log[xx])-.special.lgamma[aa];
  sm:1%aa;
  term:1%aa;
  nn:1;
  while[nn<200;
    term:term*xx%(aa+nn);
    sm:sm+term;
    if[(abs term)<1e-12*abs sm; nn:200];
    nn+:1;
  ];
  (exp lnpre)*sm
 };

/ @desc Regularized incomplete gamma Q(a,x) - Continued fraction helper
/ @param aa:float — shape parameter a > 0
/ @param xx:float — upper limit of integration x >= 0
/ @return float — Q(a,x) = 1-P(a,x) via CF
.special.gammainc_cf:{[aa;xx]
  / Continued fraction for Q(a,x) = 1-P(a,x) using Lentz's method
  / Q(a,x) = exp(-x + a*ln(x) - lgamma(a)) * CF
  lnpre:(neg xx)+(aa*log[xx])-.special.lgamma[aa];
  /
  tiny:1e-30;
  bb_cf:xx+1-aa;
  cc:1%tiny;
  dd:1%bb_cf;
  ff:dd;
  /
  nn:1;
  while[nn<=200;
    an:neg nn*(nn-aa);
    bb_cf:bb_cf+2;
    dd:bb_cf+an*dd;
    if[(abs dd)<tiny; dd:tiny];
    cc:bb_cf+an%cc;
    if[(abs cc)<tiny; cc:tiny];
    dd:1%dd;
    delta:dd*cc;
    ff:ff*delta;
    if[1e-12>abs delta-1; nn:201];
    nn+:1;
  ];
  /
  (exp lnpre)*ff
 };

/ @desc Regularized incomplete gamma function P(a,x)
/ @param aa:float — shape parameter a > 0
/ @param xx:float — upper limit of integration x >= 0
/ @return float — P(a,x)
.special.gammainc:{[aa;xx]
  if[xx<=0f; :0f];
  if[xx<aa+1;
    / Series expansion: P(a,x) = exp(-x + a*ln(x) - lgamma(a)) * sum
    :.special.gammainc_series[aa;xx]
  ];
  / Continued fraction: Q(a,x) = 1-P(a,x)
  :1-.special.gammainc_cf[aa;xx]
 };
