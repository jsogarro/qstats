/ test_nonparam_validation.q -- Signal-trap tests for Wave 5 nonparametric
/ test functions.

-1 "\n--- Nonparametric Tests Validation ---";

-1 "\n.htest.ks:";
.tst.assert_throws["ks empty x";{.htest.ks[`float$();1 2 3f]}];
.tst.assert_throws["ks empty y";{.htest.ks[1 2 3f;`float$()]}];

-1 "\n.htest.mannwhitney:";
.tst.assert_throws["mannwhitney empty x";{.htest.mannwhitney[`float$();1 2 3f]}];
.tst.assert_throws["mannwhitney empty y";{.htest.mannwhitney[1 2 3f;`float$()]}];

-1 "\n.htest.wilcoxon:";
.tst.assert_throws["wilcoxon length mismatch";{.htest.wilcoxon[1 2 3f;1 2f]}];
.tst.assert_throws["wilcoxon all zeros";{.htest.wilcoxon[1 2 3f;1 2 3f]}];

-1 "\n.htest.shapiro:";
.tst.assert_throws["shapiro n<4";{.htest.shapiro 1 2 3f}];

-1 "\n.htest.jarque_bera:";
.tst.assert_throws["jarque_bera n<4";{.htest.jarque_bera 1 2 3f}];

-1 "\n.htest nonparametric sanity checks:";

/ KS on identical samples: D=0, p=1
r1:.htest.ks[1 2 3 4 5f;1 2 3 4 5f];
.tst.assert_approx["ks identical samples D=0";r1`statistic;0f;1e-12];
.tst.assert_approx["ks identical samples p=1";r1`p_value;1f;1e-12];

/ Mann-Whitney on identical samples: U at mean, p ~ 1
r2:.htest.mannwhitney[1 2 3 4 5f;1 2 3 4 5f];
.tst.assert["mannwhitney identical samples returns dict";6=count r2];

/ Wilcoxon on identical paired samples raises (all zeros)
/ already covered above.

/ Jarque-Bera on perfectly symmetric data: JB ~ small
r3:.htest.jarque_bera (`float$1+til 100);
.tst.assert["jarque_bera symmetric data has small JB";(r3`statistic)<10f];

/ Shapiro-Wilk: just sanity-check the dict shape and W in (0, 1]
r4:.htest.shapiro `float$1+til 30;
.tst.assert["shapiro returns 6-key dict";6=count r4];
.tst.assert["shapiro W in (0, 1]";((r4`statistic)>0f) and (r4`statistic)<=1f];
