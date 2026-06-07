/ test_descriptive.q — Tie-out tests for Wave 3 descriptive statistics.
/ Compares .desc.* outputs against scipy/numpy reference values stored in
/ tests/reference/descriptive.json.

-1 "\n--- Descriptive Statistics Tie-Out Tests ---";

ref:.j.k raze read0 `:tests/reference/descriptive.json;

/ -----------------------------------------------------------------------------
/ Run each case in a list, calling `fn` with the case dict and asserting the
/ result is within `tol` of the case's `expected`. Lambdas take [c] to avoid
/ the parser ambiguity of `x`x` inside a default-param lambda.
/ -----------------------------------------------------------------------------
.tst.run_scalar_cases:{[name;cases;fn;tol]
  i:0;
  while[i<count cases;
    c:cases i;
    ay:fn c;
    .tst.assert_approx[name,"[",string[i],"]";`float$ay;`float$c`expected;tol];
    i+:1];
 };

/=============================================================================
/ 1. CENTRAL TENDENCY
/=============================================================================
-1 "\nCentral tendency:";

.tst.run_scalar_cases["mean";ref`mean;{[c] .desc.mean c`x};1e-12];
.tst.run_scalar_cases["median";ref`median;{[c] .desc.median c`x};1e-12];
.tst.run_scalar_cases["trimmed_mean";ref`trimmed_mean;{[c] .desc.trimmed_mean[c`x;c`trim]};1e-12];

/ Mode returns a list; compare sorted floats.
mode_cases:ref`mode;
{[c;i]
  ay:`float$.desc.mode c`x;
  ey:`float$c`expected;
  .tst.assert_equal["mode[",string[i],"]";ay;ey]
 }'[mode_cases;til count mode_cases];

/=============================================================================
/ 2. DISPERSION
/=============================================================================
-1 "\nDispersion:";

.tst.run_scalar_cases["var";ref`var;{[c] .desc.var c`x};1e-12];
.tst.run_scalar_cases["sd";ref`sd;{[c] .desc.sd c`x};1e-12];
.tst.run_scalar_cases["iqr";ref`iqr;{[c] .desc.iqr c`x};1e-12];
.tst.run_scalar_cases["mad";ref`mad;{[c] .desc.mad c`x};1e-12];
.tst.run_scalar_cases["cv";ref`cv;{[c] .desc.cv c`x};1e-12];

/=============================================================================
/ 3. SHAPE
/=============================================================================
-1 "\nShape:";

.tst.run_scalar_cases["skewness";ref`skewness;{[c] .desc.skewness c`x};1e-12];
.tst.run_scalar_cases["kurtosis";ref`kurtosis;{[c] .desc.kurtosis c`x};1e-12];

/=============================================================================
/ 4. ORDER STATISTICS
/=============================================================================
-1 "\nOrder statistics:";

.tst.run_scalar_cases["quantile";ref`quantile;{[c] .desc.quantile_scalar[c`x;c`p]};1e-12];
.tst.run_scalar_cases["percentile";ref`percentile;{[c] .desc.quantile_scalar[c`x;(c`p)%100f]};1e-12];

/=============================================================================
/ 5. CORRELATION
/=============================================================================
-1 "\nCorrelation:";

.tst.run_scalar_cases["cor";ref`cor;{[c] .desc.cor[c`x;c`y]};1e-12];
.tst.run_scalar_cases["spearman";ref`spearman;{[c] .desc.spearman[c`x;c`y]};1e-12];
.tst.run_scalar_cases["kendall";ref`kendall;{[c] .desc.kendall[c`x;c`y]};1e-12];

/=============================================================================
/ 6. COVARIANCE AND MATRICES
/=============================================================================
-1 "\nCovariance and matrices:";

.tst.run_scalar_cases["cov";ref`cov;{[c] .desc.cov[c`x;c`y]};1e-12];

cov_case:ref`cov_matrix;
cov_ay:.desc.cov_matrix cov_case`X;
cov_ey:cov_case`expected;
.tst.assert_approx["cov_matrix";`float$raze cov_ay;`float$raze cov_ey;1e-12];

cor_case:ref`cor_matrix;
cor_ay:.desc.cor_matrix cor_case`X;
cor_ey:cor_case`expected;
.tst.assert_approx["cor_matrix";`float$raze cor_ay;`float$raze cor_ey;1e-12];
