/ test_nonparam.q -- Tie-out tests for Wave 5 nonparametric tests.
/ Compares .htest.{ks, mannwhitney, wilcoxon, shapiro, jarque_bera}
/ against scipy reference values in tests/reference/nonparam.json.
/ Tolerances vary by test:
/  - jarque_bera: 1e-10 (closed form)
/  - mannwhitney, wilcoxon: 1e-10 statistic, 1e-6 p-value (pnorm floor)
/  - ks: 1e-10 statistic, 5e-6 p-value (asymptotic series truncation)
/  - shapiro: 1e-4 statistic, 1e-2 p-value (Royston 1992 approximation)

-1 "\n--- Nonparametric Tests Tie-Out ---";

ref:.j.k raze read0 `:tests/reference/nonparam.json;

/ ----- KS -----
-1 "\nKolmogorov-Smirnov:";
ks_cases:ref`ks;
{[c;i]
  res:.htest.ks[c`x;c`y];
  .tst.assert_approx["ks[",string[i],"] D";`float$res`statistic;`float$c`statistic;1e-10];
  .tst.assert_approx["ks[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;5e-6]
 }'[ks_cases;til count ks_cases];

/ ----- Mann-Whitney -----
-1 "\nMann-Whitney U:";
mw_cases:ref`mannwhitney;
{[c;i]
  res:.htest.mannwhitney[c`x;c`y];
  .tst.assert_approx["mannwhitney[",string[i],"] U";`float$res`statistic;`float$c`statistic;1e-10];
  .tst.assert_approx["mannwhitney[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;1e-6]
 }'[mw_cases;til count mw_cases];

/ ----- Wilcoxon signed-rank -----
/ scipy returns min(W+, W-) for two-sided; our impl returns W+. Compare
/ p-values primarily; for the statistic compare against w_plus.
-1 "\nWilcoxon signed-rank:";
wcx_cases:ref`wilcoxon;
{[c;i]
  res:.htest.wilcoxon[c`x;c`y];
  .tst.assert_approx["wilcoxon[",string[i],"] W+";`float$res`statistic;`float$c`w_plus;1e-10];
  .tst.assert_approx["wilcoxon[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;1e-6]
 }'[wcx_cases;til count wcx_cases];

/ ----- Shapiro-Wilk -----
/ Royston 1992 approximation; tolerance loosened. n < 12 returns p=0n.
-1 "\nShapiro-Wilk:";
sw_cases:ref`shapiro;
{[c;i]
  res:.htest.shapiro c`x;
  .tst.assert_approx["shapiro[",string[i],"] W";`float$res`statistic;`float$c`statistic;1e-4];
  if[(count c`x)>=12;
    .tst.assert_approx["shapiro[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;0.01]]
 }'[sw_cases;til count sw_cases];

/ ----- Shapiro-Wilk small-sample (n in [4, 11]) -----
/ Wave 8 Gap B2: Royston small-sample transform
-1 "\nShapiro-Wilk small-sample (n in [4, 11]):";
sw_small_cases:ref`shapiro_small;
{[c;i]
  res:.htest.shapiro c`x;
  .tst.assert_approx["shapiro_small[",string[i],"] W";`float$res`statistic;`float$c`statistic;1e-4];
  / p-value tolerance for small n is documented as ~1e-2 (Royston's accuracy floor)
  / But scipy's implementation is tighter; use 5e-3 to allow for the polynomial approx.
  .tst.assert_approx["shapiro_small[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;5e-3]
 }'[sw_small_cases;til count sw_small_cases];

/ ----- Jarque-Bera -----
-1 "\nJarque-Bera:";
jb_cases:ref`jarque_bera;
{[c;i]
  res:.htest.jarque_bera c`x;
  .tst.assert_approx["jarque_bera[",string[i],"] JB";`float$res`statistic;`float$c`statistic;1e-10];
  .tst.assert_approx["jarque_bera[",string[i],"] p_value";`float$res`p_value;`float$c`p_value;1e-10]
 }'[jb_cases;til count jb_cases];
