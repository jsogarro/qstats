/=============================================================================
/ descriptive.q — Descriptive statistics
/ Namespace: .desc
/ Wave 3 of qstats master plan.
/ All functions match scipy/numpy reference within 1e-10 unless noted.
/ Skewness and kurtosis use type-2 (R/SAS/SPSS, scipy `bias=False`).
/ Quantile uses type-7 (R/numpy default; linear interpolation).
/=============================================================================

.desc.validate:{[check;msg]
  if[not check;'`$"invalid_arg: ",msg];
 };

/=============================================================================
/ 4. ORDER STATISTICS (defined first; iqr/mad/median depend on quantile)
/=============================================================================

/ @desc Type-7 quantile (R/numpy default): linear interpolation on order stats.
/ @param x:float[] — input vector (count >= 1)
/ @param p:float — probability in [0, 1]
.desc.quantile_scalar:{[x;p]
  .desc.validate[(p>=0f) and p<=1f;"p must be in [0, 1]"];
  n:count x;
  .desc.validate[n>0;"x must be non-empty"];
  s:asc x;
  if[n=1; :s 0];
  if[p>=1f; :s n-1];
  if[p<=0f; :s 0];
  h:(n-1)*p;
  hf:`long$floor h;
  fr:h-hf;
  s[hf]+fr*s[hf+1]-s[hf]
 };
.desc.quantile:{[x;p] .desc.quantile_scalar[x;] each p};

/ @desc Percentile — same as quantile but with percentage input (0-100)
.desc.percentile:{[x;p] .desc.quantile[x;p%100f]};

/=============================================================================
/ 1. CENTRAL TENDENCY
/=============================================================================

.desc.mean:{[x] avg x};

.desc.median:{[x] .desc.quantile_scalar[x;0.5]};

/ Mode: most frequent value(s); may be multi-modal. Returns ascending list.
/ `where dict` returns the dict's keys whose values are truthy, so no need
/ to apply `key` separately (doing so double-indexes and trips a type error).
.desc.mode:{[x]
  fr:count each group x;
  asc where fr=max fr
 };

/ Trimmed mean: discard `tr` fraction from each tail, then average.
/ Parameter is named `tr` (not `trim`) because `trim` is a q builtin.
.desc.trimmed_mean:{[x;tr]
  .desc.validate[(tr>=0f) and tr<0.5;"trim fraction must be in [0, 0.5)"];
  n:count x;
  k:`long$floor tr*n;
  avg (asc x)[k+til n-2*k]
 };

/=============================================================================
/ 2. DISPERSION
/=============================================================================

/ Sample variance (unbiased, n-1 denominator).
.desc.var:{[x]
  n:count x;
  .desc.validate[n>1;"variance requires at least 2 observations"];
  m:avg x;
  (sum (x-m) xexp 2)%n-1
 };

.desc.sd:{[x] sqrt .desc.var x};

.desc.iqr:{[x] (.desc.quantile_scalar[x;0.75])-.desc.quantile_scalar[x;0.25]};

/ Median absolute deviation, scaled to be a consistent SD estimator for
/ Gaussian data. Scaling constant is 1 / qnorm(0.75) to full double precision
/ to match scipy.stats.median_abs_deviation(scale='normal'). R's `mad` uses
/ the rounded 1.4826 by default — that floor is too coarse for our 1e-12
/ tieout tolerance.
.desc.mad:{[x]
  m:.desc.quantile_scalar[x;0.5];
  1.4826022185056018*.desc.quantile_scalar[abs x-m;0.5]
 };

/ Coefficient of variation: sd / |mean|. Useful for cross-scale comparison.
.desc.cv:{[x]
  m:avg x;
  .desc.validate[not m=0f;"cv undefined when mean is zero"];
  (.desc.sd x)%abs m
 };

/=============================================================================
/ 3. SHAPE
/=============================================================================

/ Sample skewness, type-2 (R/SAS/SPSS default; scipy `bias=False`).
/ G1 = sqrt(n*(n-1))/(n-2) * (m3 / m2^(3/2)) with 1/n moments.
.desc.skewness:{[x]
  n:count x;
  .desc.validate[n>=3;"skewness requires at least 3 observations"];
  m:avg x;
  m2:(sum (x-m) xexp 2)%n;
  m3:(sum (x-m) xexp 3)%n;
  g1:m3%m2 xexp 1.5;
  g1*(sqrt n*n-1f)%n-2f
 };

/ Excess kurtosis, type-2 (R/SAS/SPSS; scipy `bias=False, fisher=True`).
/ G2 = ((n-1)/((n-2)(n-3))) * ((n+1)*g2 + 6), g2 = m4/m2^2 - 3.
.desc.kurtosis:{[x]
  n:count x;
  .desc.validate[n>=4;"kurtosis requires at least 4 observations"];
  m:avg x;
  m2:(sum (x-m) xexp 2)%n;
  m4:(sum (x-m) xexp 4)%n;
  g2:(m4%m2 xexp 2)-3f;
  ((n-1f)%(n-2f)*n-3f)*((n+1f)*g2)+6f
 };

/=============================================================================
/ 5. CORRELATION
/=============================================================================

.desc.cor:{[x;y]
  .desc.validate[(count x)=count y;"x and y must have the same length"];
  mx:avg x;
  my:avg y;
  dx:x-mx;
  dy:y-my;
  (sum dx*dy)%(sqrt sum dx xexp 2)*sqrt sum dy xexp 2
 };

/ Spearman rank correlation (Pearson on ranks; handles ties via average rank).
.desc.spearman:{[x;y]
  .desc.validate[(count x)=count y;"x and y must have the same length"];
  rx:`float$iasc iasc x;
  ry:`float$iasc iasc y;
  .desc.cor[rx;ry]
 };

/ Kendall tau-b (handles ties). O(n^2) reference implementation.
.desc.kendall:{[x;y]
  .desc.validate[(count x)=count y;"x and y must have the same length"];
  n:count x;
  / Pairwise sign matrices. sx[i;j] = sign(x[i] - x[j]).
  sx:signum (x-/:x);
  sy:signum (y-/:y);
  / Sum of sx*sy over i<j pairs. Full matrix sums each pair twice (i,j and
  / j,i) with sign(sx*sy) symmetric, so halve to get upper-triangle count.
  num:0.5*sum sum sx*sy;
  / Tie counts: number of unordered pairs (i,j) with i!=j where x[i]=x[j].
  tx:0.5*sum sum (sx=0f) and not (til n)=/:til n;
  ty:0.5*sum sum (sy=0f) and not (til n)=/:til n;
  n_pairs:0.5*n*n-1;
  num%sqrt (n_pairs-tx)*n_pairs-ty
 };

/=============================================================================
/ 6. COVARIANCE AND CORRELATION MATRICES
/=============================================================================

.desc.cov:{[x;y]
  .desc.validate[(count x)=count y;"x and y must have the same length"];
  n:count x;
  .desc.validate[n>1;"covariance requires at least 2 observations"];
  mx:avg x;
  my:avg y;
  (sum (x-mx)*y-my)%n-1
 };

/ Covariance matrix of a data matrix X (n rows by p columns).
.desc.cov_matrix:{[X]
  n:count X;
  .desc.validate[n>1;"cov_matrix requires at least 2 rows"];
  means:avg each flip X;
  Xc:X-\:means;
  (.la.crossprod Xc)%n-1
 };

/ Pearson correlation matrix derived from the covariance matrix.
.desc.cor_matrix:{[X]
  C:.desc.cov_matrix X;
  sds:sqrt .la.diag C;
  C%sds*/:\:sds
 };

/=============================================================================
/ 7. SUMMARY STATISTICS
/=============================================================================

/ 8-number summary dictionary: n, mean, sd, min, q25, median, q75, max.
.desc.summary:{[x]
  `n`mean`sd`min`q25`median`q75`max!(
    count x;
    avg x;
    .desc.sd x;
    min x;
    .desc.quantile_scalar[x;0.25];
    .desc.quantile_scalar[x;0.5];
    .desc.quantile_scalar[x;0.75];
    max x)
 };

/ Frequency table: val | cnt | prop. Column is `val` rather than `value`
/ because `value` is a q reserved word (the dict-value builtin), which q's
/ table literal syntax cannot bind to as a column name.
.desc.freq:{[x]
  g:count each group x;
  vals:asc key g;
  cnts:g vals;
  tot:sum cnts;
  ([] val:vals; cnt:cnts; prop:cnts%tot)
 };

/ Cross-tabulation (contingency table) of two categorical vectors.
/ Returns a matrix indexed by (asc distinct x) x (asc distinct y).
.desc.crosstab:{[x;y]
  .desc.validate[(count x)=count y;"x and y must have the same length"];
  xv:asc distinct x;
  yv:asc distinct y;
  g:count each group flip (x;y);
  cell:{[g;xi;yj] $[(xi;yj) in key g; g (xi;yj); 0]};
  {[cell;g;yv;xi] cell[g;xi;] each yv}[cell;g;yv;] each xv
 };
