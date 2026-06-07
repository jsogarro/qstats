/=============================================================================
/ tests.q — Parametric hypothesis tests
/ Namespace: .htest (not .test, to avoid collision with the test framework)
/ Wave 4 of qstats master plan.
/ Every test returns a dictionary with these keys:
/   statistic   — test statistic
/   df          — degrees of freedom (atom, pair, or null)
/   p_value     — two-sided p-value
/   method      — string description
/   alternative — "two.sided" | "less" | "greater"
/   ci          — 2-element list (lo, hi); (0n;0n) when not computed
/=============================================================================

.htest.validate:{[check;msg]
  if[not check;'`$"invalid_arg: ",msg];
 };

/=============================================================================
/ 1. t-Tests
/=============================================================================

/ One-sample t-test: t = (mx - mu0) / (s / sqrt(n)), df = n - 1.
.htest.ttest1:{[x;mu0]
  n:count x;
  .htest.validate[n>1;"ttest1 requires at least 2 observations"];
  mx:avg x;
  s:.desc.sd x;
  se:s%sqrt n;
  tstat:(mx-mu0)%se;
  df:n-1;
  pval:2f*.dist.pt[neg abs tstat;df];
  tcrit:.dist.qt[0.975;df];
  ci:(mx-tcrit*se;mx+tcrit*se);
  `statistic`df`p_value`method`alternative`ci!(tstat;df;pval;"One Sample t-test";"two.sided";ci)
 };

/ Two-sample t-test, equal variance (pooled SD), df = nx + ny - 2.
.htest.ttest2:{[x;y]
  nx:count x;
  ny:count y;
  .htest.validate[(nx>1) and ny>1;"ttest2 requires at least 2 observations per group"];
  mx:avg x;
  my:avg y;
  vx:.desc.var x;
  vy:.desc.var y;
  sp:sqrt (((nx-1)*vx)+(ny-1)*vy)%nx+ny-2;
  se:sp*sqrt (1f%nx)+1f%ny;
  tstat:(mx-my)%se;
  df:nx+ny-2;
  pval:2f*.dist.pt[neg abs tstat;df];
  `statistic`df`p_value`method`alternative`ci!(tstat;df;pval;"Two Sample t-test (equal variance)";"two.sided";(0n;0n))
 };

/ Welch's t-test (unequal variance), Welch-Satterthwaite df.
.htest.welch:{[x;y]
  nx:count x;
  ny:count y;
  .htest.validate[(nx>1) and ny>1;"welch requires at least 2 observations per group"];
  mx:avg x;
  my:avg y;
  vx:.desc.var x;
  vy:.desc.var y;
  vx_n:vx%nx;
  vy_n:vy%ny;
  se:sqrt vx_n+vy_n;
  tstat:(mx-my)%se;
  num:(vx_n+vy_n) xexp 2;
  den:((vx_n xexp 2)%nx-1)+(vy_n xexp 2)%ny-1;
  df:num%den;
  pval:2f*.dist.pt[neg abs tstat;df];
  `statistic`df`p_value`method`alternative`ci!(tstat;df;pval;"Welch Two Sample t-test";"two.sided";(0n;0n))
 };

/ Paired t-test: reduces to one-sample t-test on differences.
.htest.ttest_paired:{[x;y]
  .htest.validate[(count x)=count y;"x and y must have the same length"];
  res:.htest.ttest1[x-y;0f];
  res[`method]:"Paired t-test";
  res
 };

/=============================================================================
/ 2. F-Test for Variance Ratio
/=============================================================================

/ F-test for equality of variances: F = var(x) / var(y).
/ Two-sided p-value via the standard min(cdf, 1-cdf) construction.
.htest.ftest:{[x;y]
  nx:count x;
  ny:count y;
  .htest.validate[(nx>1) and ny>1;"ftest requires at least 2 observations per group"];
  vx:.desc.var x;
  vy:.desc.var y;
  fstat:vx%vy;
  df1:nx-1;
  df2:ny-1;
  cdf:.dist.pf[fstat;df1;df2];
  pval:2f*$[cdf<0.5;cdf;1f-cdf];
  `statistic`df`p_value`method`alternative`ci!(fstat;(df1;df2);pval;"F test to compare two variances";"two.sided";(0n;0n))
 };

/=============================================================================
/ 3. Chi-Squared Tests
/=============================================================================

/ Chi-squared goodness-of-fit: sum (O - E)^2 / E, df = k - 1.
.htest.chisq_gof:{[observed;expected]
  .htest.validate[(count observed)=count expected;"observed and expected must have the same length"];
  k:count observed;
  .htest.validate[k>1;"chisq_gof requires at least 2 categories"];
  chi:sum ((observed-expected) xexp 2)%expected;
  df:k-1;
  pval:1f-.dist.pchisq[chi;df];
  `statistic`df`p_value`method`alternative`ci!(chi;df;pval;"Chi-squared goodness-of-fit test";"two.sided";(0n;0n))
 };

/ Chi-squared test of independence on a 2-D contingency table.
/ E_ij = (row_i_total * col_j_total) / grand_total. df = (r-1)(c-1).
/ For 2x2 tables, Yates' continuity correction is applied (matching R's
/ chisq.test and scipy's chi2_contingency defaults): replace |O - E| with
/ max(0, |O - E| - 0.5) before squaring.
.htest.chisq_ind:{[tbl]
  r:count tbl;
  c:count first tbl;
  .htest.validate[(r>1) and c>1;"chisq_ind requires at least a 2x2 table"];
  rowtot:sum each tbl;
  coltot:sum tbl;
  grand:sum rowtot;
  expected:(rowtot*/:\:coltot)%grand;
  diff:abs tbl-expected;
  if[(r=2) and c=2; diff:0f|diff-0.5];
  chi:sum raze (diff xexp 2)%expected;
  df:(r-1)*c-1;
  pval:1f-.dist.pchisq[chi;df];
  `statistic`df`p_value`method`alternative`ci!(chi;df;pval;"Pearson's Chi-squared test of independence";"two.sided";(0n;0n))
 };

/=============================================================================
/ 4. One-Way ANOVA
/=============================================================================

/ One-way ANOVA F-test. y is the response vector; grp is the parallel group
/ vector (symbols, longs, anything hashable). F = MS_between / MS_within.
.htest.anova1:{[y;grp]
  .htest.validate[(count y)=count grp;"y and grp must have the same length"];
  n:count y;
  k:count distinct grp;
  .htest.validate[k>1;"anova1 requires at least 2 groups"];
  .htest.validate[n>k;"anova1 requires more observations than groups"];
  my:avg y;
  g:group grp;
  grp_means:value avg each y g;
  grp_cnts:value count each g;
  ssbetween:sum grp_cnts*(grp_means-my) xexp 2;
  sstotal:sum (y-my) xexp 2;
  sswithin:sstotal-ssbetween;
  df1:k-1;
  df2:n-k;
  fstat:(ssbetween%df1)%sswithin%df2;
  pval:1f-.dist.pf[fstat;df1;df2];
  `statistic`df`p_value`method`alternative`ci!(fstat;(df1;df2);pval;"One-way ANOVA";"two.sided";(0n;0n))
 };

/=============================================================================
/ 5. Correlation Test (Pearson)
/=============================================================================

/ t = r * sqrt(n - 2) / sqrt(1 - r^2), df = n - 2.
.htest.cortest:{[x;y]
  .htest.validate[(count x)=count y;"x and y must have the same length"];
  n:count x;
  .htest.validate[n>2;"cortest requires at least 3 observations"];
  r:.desc.cor[x;y];
  tstat:(r*sqrt n-2)%sqrt 1f-r*r;
  df:n-2;
  pval:2f*.dist.pt[neg abs tstat;df];
  `statistic`df`p_value`method`alternative`ci!(tstat;df;pval;"Pearson's product-moment correlation";"two.sided";(0n;0n))
 };

/=============================================================================
/ Wave 5 helpers
/=============================================================================

/ Tie-averaged ranks (1-indexed). For ties, all tied positions get the mean
/ of the would-be sequential ranks. Matches scipy.stats.rankdata default.
.htest._ranks:{[xs]
  raw:1f+iasc iasc xs;
  g:group xs;
  mean_per_val:{[r;ix] avg r ix}[raw;] each g;
  mean_per_val xs
 };

/=============================================================================
/ 6. One-Proportion Test
/=============================================================================

/ One-sample z-test for a proportion (normal approximation, no continuity
/ correction). z = (phat - p0) / sqrt(p0 (1 - p0) / n).
.htest.proptest:{[xx;n;p0]
  .htest.validate[(p0>0f) and p0<1f;"p0 must be in (0, 1)"];
  .htest.validate[n>0;"n must be positive"];
  .htest.validate[(xx>=0) and xx<=n;"successes must be in [0, n]"];
  phat:xx%n;
  se:sqrt (p0*1f-p0)%n;
  zstat:(phat-p0)%se;
  pval:2f*.dist.pnorm[neg abs zstat;0f;1f];
  `statistic`df`p_value`method`alternative`ci!(zstat;0n;pval;"1-sample test for proportion (no continuity correction)";"two.sided";(0n;0n))
 };

/=============================================================================
/ 7. Kolmogorov-Smirnov Test (Wave 5)
/=============================================================================

/ Two-sample KS test. D = max |F_x(t) - F_y(t)| over the pooled support.
/ Asymptotic p-value via Kolmogorov distribution series with Stephens (1970)
/ finite-sample correction lambda = (en + 0.12 + 0.11/en) * D, which is what
/ scipy.stats.ks_2samp(method='asymp') uses for the p-value.
.htest.ks:{[x;y]
  nx:count x;
  ny:count y;
  .htest.validate[(nx>0) and ny>0;"both samples must be non-empty"];
  pts:asc distinct x,y;
  fx:{[xs;t] (sum xs<=t)%count xs}[x;] each pts;
  fy:{[ys;t] (sum ys<=t)%count ys}[y;] each pts;
  d:max abs fx-fy;
  en:sqrt (nx*ny)%nx+ny;
  / scipy.ks_2samp(method='asymp') uses lam = en * D (no Stephens correction).
  / Edge case: D=0 (samples have identical empirical CDFs) -> SF(0) = 1.
  / The asymptotic series sum (-1)^(k-1) at lam=0 doesn't converge, so handle
  / explicitly.
  lam:en*d;
  pval:$[d=0f;
    1f;
    [ks_terms:1+til 100;
     signs:`float$xexp[-1f;ks_terms-1];
     0f|1f&2f*sum signs*exp neg 2f*(ks_terms*ks_terms)*lam*lam]];
  `statistic`df`p_value`method`alternative`ci!(d;0n;pval;"Two-sample Kolmogorov-Smirnov test";"two.sided";(0n;0n))
 };

/=============================================================================
/ 8. Mann-Whitney U Test (Wave 5)
/=============================================================================

/ Mann-Whitney U (Wilcoxon rank-sum). Two-sided normal approximation with
/ tie correction and continuity correction, matching scipy.stats.mannwhitneyu
/ with method='asymptotic' and use_continuity=True. Returns U1 (sum of ranks
/ of x minus the offset), which is what scipy reports.
.htest.mannwhitney:{[x;y]
  nx:count x;
  ny:count y;
  .htest.validate[(nx>0) and ny>0;"both samples must be non-empty"];
  Nt:nx+ny;
  combined:x,y;
  ranks:.htest._ranks combined;
  Rx:sum nx#ranks;
  U1:Rx-(nx*nx+1f)%2f;
  mu:(nx*ny)%2f;
  / Variance with tie correction
  g:group combined;
  tsz:`float$value count each g;
  tie_adj:sum (tsz*tsz*tsz)-tsz;
  var_u:((nx*ny)%12f)*((Nt+1f)-tie_adj%Nt*Nt-1f);
  / Continuity correction
  z:((abs U1-mu)-0.5)%sqrt var_u;
  pval:2f*1f-.dist.pnorm[abs z;0f;1f];
  `statistic`df`p_value`method`alternative`ci!(U1;0n;pval;"Mann-Whitney U test (Wilcoxon rank-sum)";"two.sided";(0n;0n))
 };

/=============================================================================
/ 9. Wilcoxon Signed-Rank Test (Wave 5)
/=============================================================================

/ Wilcoxon signed-rank for paired samples. Drops zero differences, ranks
/ |d| with ties averaged, computes W+. Two-sided normal approximation with
/ tie correction and continuity correction, matching
/ scipy.stats.wilcoxon(method='approx', correction=True). Returns W+ as the
/ statistic (scipy reports min(W+, W-); we adjust in the tieout test).
.htest.wilcoxon:{[x;y]
  .htest.validate[(count x)=count y;"x and y must have the same length"];
  d:x-y;
  d:d where not d=0f;
  n:count d;
  .htest.validate[n>0;"all paired differences are zero"];
  ad:abs d;
  ranks:.htest._ranks ad;
  W_plus:sum ranks where d>0f;
  W_minus:sum ranks where d<0f;
  mu:(n*n+1f)%4f;
  var_w:(n*(n+1f)*(2f*n)+1f)%24f;
  / Tie correction: subtract sum(t^3 - t)/48
  g:group ad;
  tsz:`float$value count each g;
  tie_adj:(sum (tsz*tsz*tsz)-tsz)%48f;
  var_w:var_w-tie_adj;
  / Continuity correction
  z:((abs W_plus-mu)-0.5)%sqrt var_w;
  pval:2f*1f-.dist.pnorm[abs z;0f;1f];
  `statistic`df`p_value`method`alternative`ci!(W_plus;0n;pval;"Wilcoxon signed-rank test";"two.sided";(0n;0n))
 };

/=============================================================================
/ 10. Jarque-Bera Normality Test (Wave 5)
/=============================================================================

/ JB = n/6 * (S^2 + (K - 3)^2 / 4) using BIASED (type-1) skewness and kurtosis
/ moments, matching scipy.stats.jarque_bera. Under H0, JB ~ chi^2(2).
.htest.jarque_bera:{[xx]
  n:count xx;
  .htest.validate[n>=4;"jarque_bera requires at least 4 observations"];
  m:avg xx;
  d:xx-m;
  m2:(sum d xexp 2)%n;
  m3:(sum d xexp 3)%n;
  m4:(sum d xexp 4)%n;
  S:m3%m2 xexp 1.5;
  K_excess:(m4%m2 xexp 2)-3f;
  jb:(n%6f)*(S*S)+(K_excess*K_excess)%4f;
  pval:1f-.dist.pchisq[jb;2];
  `statistic`df`p_value`method`alternative`ci!(jb;2;pval;"Jarque-Bera normality test";"two.sided";(0n;0n))
 };

/=============================================================================
/ 11. Shapiro-Wilk Normality Test (Wave 5)
/=============================================================================

/ Shapiro-Wilk W via Royston (1992) weight construction (Algorithm AS R94).
/ For n in [12, 5000] the p-value uses Royston's log-normal transformation
/ of log(1 - W); see Royston 1992 "Approximating the Shapiro-Wilk W-test for
/ non-normality." Tolerance vs scipy is ~1e-4 for W, ~1e-3 for p-value at
/ moderate n; documented as a known precision floor.
.htest.shapiro:{[xx]
  n:count xx;
  .htest.validate[(n>=4) and n<=5000;"shapiro requires 4 <= n <= 5000"];
  sx:asc xx;
  mn:avg xx;
  s2:sum (xx-mn) xexp 2;
  / m_i = Phi^-1((i - 3/8) / (n + 1/4))
  i_vec:1f+til n;
  m:.dist.qnorm[(i_vec-0.375)%n+0.25;0f;1f];
  cc:sqrt sum m*m;
  uu:1f%sqrt n;
  / Royston polynomial coefficients for a_n and a_{n-1}
  a_n:(((((neg 2.706056)*uu xexp 5)+4.434685*uu xexp 4)-2.071190*uu xexp 3)-0.147981*uu xexp 2)+(0.221157*uu)+m[n-1]%cc;
  a_n1:(((((neg 3.582633)*uu xexp 5)+5.682633*uu xexp 4)-1.752460*uu xexp 3)-0.293762*uu xexp 2)+(0.042981*uu)+m[n-2]%cc;
  / Solve for e^2 using m^T m - 2 m_n^2 - 2 m_{n-1}^2 and 1 - 2 a_n^2 - 2 a_{n-1}^2
  mm:sum m*m;
  num:(mm-2f*m[n-1]*m[n-1])-2f*m[n-2]*m[n-2];
  den:(1f-2f*a_n*a_n)-2f*a_n1*a_n1;
  e_:sqrt num%den;
  / Build a-vector: a[0] = -a_n, a[1] = -a_{n-1}, a[n-2] = a_{n-1}, a[n-1] = a_n,
  / a[i] = m_i / e for i in 2..n-3 (0-indexed).
  a:n#0f;
  a[0]:neg a_n;
  a[1]:neg a_n1;
  a[n-2]:a_n1;
  a[n-1]:a_n;
  mid_idx:2+til n-4;
  a[mid_idx]:m[mid_idx]%e_;
  / W statistic
  num_w:sum a*sx;
  W:(num_w*num_w)%s2;
  / p-value via Royston transformation. For n in [12, 5000] use log(1-W) form;
  / for n in [4, 11] use a polynomial form (omitted here -- accuracy lower).
  / Royston 1992 coefficients (n in [12, 5000]):
  /   mu  = 0.0038915*g^3 - 0.083751*g^2 - 0.31082*g - 1.5861
  /   sig = exp(0.0030302*g^2 - 0.082676*g - 0.4803)
  / Parenthesise each subtraction left-to-right so q's right-associative `-`
  / doesn't flip signs (the trap that initially gave us mu = ... + 1.5861).
  pval:$[n>=12;
    [g:log n;
     mu_y:(((0.0038915*g xexp 3)-0.083751*g*g)-0.31082*g)-1.5861;
     sig_y:exp ((0.0030302*g*g)-0.082676*g)-0.4803;
     yval:log 1f-W;
     z:(yval-mu_y)%sig_y;
     1f-.dist.pnorm[z;0f;1f]];
    / n < 12 returns 0n; small-sample p-value form not implemented
    0n];
  `statistic`df`p_value`method`alternative`ci!(W;0n;pval;"Shapiro-Wilk normality test";"two.sided";(0n;0n))
 };
