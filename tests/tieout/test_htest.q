/ test_htest.q -- Tie-out tests for Wave 4 parametric hypothesis tests.
/ Compares .htest.* outputs against scipy reference values in
/ tests/reference/htest.json. Tests are structured as:
/   { inputs..., statistic, p_value [, df, df1, df2, r] }
/ For each case we verify the test statistic and p-value match within tol,
/ and (where applicable) the degrees of freedom.

-1 "\n--- Parametric Hypothesis Tests Tie-Out ---";

ref:.j.k raze read0 `:tests/reference/htest.json;

/ ---------------------------------------------------------------------------
/ Helper: run a list of cases, applying `fn` to each case dict and asserting
/ the returned result dict matches the case's `statistic` and `p_value` keys
/ within `tol` (and `df` if `check_df` is 1b).
/ ---------------------------------------------------------------------------
.tst.run_htest_cases:{[name;cases;fn;tol;check_df]
  i:0;
  while[i<count cases;
    c:cases i;
    res:fn c;
    .tst.assert_approx[name,"[",string[i],"] statistic";`float$res`statistic;`float$c`statistic;tol];
    .tst.assert_approx[name,"[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;tol];
    if[check_df;
      .tst.assert_approx[name,"[",string[i],"] df";`float$res`df;`float$c`df;1e-10]];
    i+:1];
 };

/=============================================================================
/ 1. t-Tests
/=============================================================================
-1 "\nt-tests:";

.tst.run_htest_cases["ttest1";ref`ttest1;{[c] .htest.ttest1[c`x;c`mu0]};1e-10;1b];

/ ttest2 CI
{[c;i]
  res:.htest.ttest2[c`x;c`y];
  .tst.assert_approx["ttest2[",string[i],"] statistic";`float$res`statistic;`float$c`statistic;1e-10];
  .tst.assert_approx["ttest2[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;1e-10];
  .tst.assert_approx["ttest2[",string[i],"] df";`float$res`df;`float$c`df;1e-10];
  .tst.assert_approx["ttest2[",string[i],"] ci_lo";res[`ci][0];c`ci_lo;1e-9];
  .tst.assert_approx["ttest2[",string[i],"] ci_hi";res[`ci][1];c`ci_hi;1e-9]
 }'[ref`ttest2;til count ref`ttest2];

/ welch CI
{[c;i]
  res:.htest.welch[c`x;c`y];
  .tst.assert_approx["welch[",string[i],"] statistic";`float$res`statistic;`float$c`statistic;1e-10];
  .tst.assert_approx["welch[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;1e-10];
  .tst.assert_approx["welch[",string[i],"] df";`float$res`df;`float$c`df;1e-10];
  .tst.assert_approx["welch[",string[i],"] ci_lo";res[`ci][0];c`ci_lo;1e-9];
  .tst.assert_approx["welch[",string[i],"] ci_hi";res[`ci][1];c`ci_hi;1e-9]
 }'[ref`welch;til count ref`welch];

/ ttest_paired CI
{[c;i]
  res:.htest.ttest_paired[c`x;c`y];
  .tst.assert_approx["ttest_paired[",string[i],"] statistic";`float$res`statistic;`float$c`statistic;1e-10];
  .tst.assert_approx["ttest_paired[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;1e-10];
  .tst.assert_approx["ttest_paired[",string[i],"] df";`float$res`df;`float$c`df;1e-10];
  .tst.assert_approx["ttest_paired[",string[i],"] ci_lo";res[`ci][0];c`ci_lo;1e-9];
  .tst.assert_approx["ttest_paired[",string[i],"] ci_hi";res[`ci][1];c`ci_hi;1e-9]
 }'[ref`ttest_paired;til count ref`ttest_paired];

/=============================================================================
/ 2. F-Test
/=============================================================================
-1 "\nF-test:";

ftest_cases:ref`ftest;
{[c;i]
  res:.htest.ftest[c`x;c`y];
  .tst.assert_approx["ftest[",string[i],"] statistic";`float$res`statistic;`float$c`statistic;1e-10];
  .tst.assert_approx["ftest[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;1e-10];
  / df is a (df1; df2) pair; cast both sides to float so 49 (long) and 49f match.
  .tst.assert_equal["ftest[",string[i],"] df";`float$res`df;`float$(c`df1;c`df2)];
  .tst.assert_approx["ftest[",string[i],"] ci_lo";res[`ci][0];c`ci_lo;1e-9];
  .tst.assert_approx["ftest[",string[i],"] ci_hi";res[`ci][1];c`ci_hi;1e-9]
 }'[ftest_cases;til count ftest_cases];

/=============================================================================
/ 3. Chi-Squared
/=============================================================================
-1 "\nChi-squared:";

.tst.run_htest_cases["chisq_gof";ref`chisq_gof;{[c] .htest.chisq_gof[c`observed;c`expected]};1e-10;1b];

ind_cases:ref`chisq_ind;
{[c;i]
  res:.htest.chisq_ind c`table;
  .tst.assert_approx["chisq_ind[",string[i],"] statistic";`float$res`statistic;`float$c`statistic;1e-10];
  .tst.assert_approx["chisq_ind[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;1e-10];
  .tst.assert_approx["chisq_ind[",string[i],"] df";`float$res`df;`float$c`df;1e-10]
 }'[ind_cases;til count ind_cases];

/=============================================================================
/ 4. ANOVA
/=============================================================================
-1 "\nANOVA:";

anova_cases:ref`anova1;
{[c;i]
  res:.htest.anova1[c`y;c`grp];
  .tst.assert_approx["anova1[",string[i],"] F";`float$res`statistic;`float$c`statistic;1e-10];
  .tst.assert_approx["anova1[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;1e-10];
  .tst.assert_equal["anova1[",string[i],"] df";`float$res`df;`float$(c`df1;c`df2)]
 }'[anova_cases;til count anova_cases];

/=============================================================================
/ 5. Correlation test
/=============================================================================
-1 "\nCorrelation test:";

/ cortest CI
{[c;i]
  res:.htest.cortest[c`x;c`y];
  .tst.assert_approx["cortest[",string[i],"] statistic";`float$res`statistic;`float$c`statistic;1e-10];
  .tst.assert_approx["cortest[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;1e-10];
  .tst.assert_approx["cortest[",string[i],"] df";`float$res`df;`float$c`df;1e-10];
  .tst.assert_approx["cortest[",string[i],"] ci_lo";res[`ci][0];c`ci_lo;1e-9];
  .tst.assert_approx["cortest[",string[i],"] ci_hi";res[`ci][1];c`ci_hi;1e-9]
 }'[ref`cortest;til count ref`cortest];

/=============================================================================
/ 6. Proportion test
/=============================================================================
-1 "\nProportion test:";

/ proptest p-values use pnorm, which has a 5e-7 precision floor per the
/ documented Numerical Accuracy section. Loosen the tolerance accordingly.
/ proptest CI
{[c;i]
  res:.htest.proptest[c`x;c`n;c`p0];
  .tst.assert_approx["proptest[",string[i],"] statistic";`float$res`statistic;`float$c`statistic;1e-6];
  .tst.assert_approx["proptest[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;1e-6];
  .tst.assert_approx["proptest[",string[i],"] ci_lo";res[`ci][0];c`ci_lo;1e-9];
  .tst.assert_approx["proptest[",string[i],"] ci_hi";res[`ci][1];c`ci_hi;1e-9]
 }'[ref`proptest;til count ref`proptest];
