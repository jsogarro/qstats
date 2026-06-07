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
