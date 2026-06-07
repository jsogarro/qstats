/ test_htest_validation.q -- Input validation tests for .htest.* functions.
/ Tieout suite covers numeric correctness; these tests cover signal handling
/ for invalid arguments.

-1 "\n--- Parametric Tests Validation ---";

-1 "\n.htest.ttest1:";
.tst.assert_throws["ttest1 with n=1";{.htest.ttest1[enlist 5f;0f]}];

-1 "\n.htest.ttest2 / welch:";
.tst.assert_throws["ttest2 with empty x";{.htest.ttest2[`float$();1 2 3f]}];
.tst.assert_throws["ttest2 with n=1 group";{.htest.ttest2[enlist 1f;1 2 3f]}];
.tst.assert_throws["welch with n=1 group";{.htest.welch[1 2 3f;enlist 1f]}];

-1 "\n.htest.ttest_paired:";
.tst.assert_throws["ttest_paired length mismatch";{.htest.ttest_paired[1 2 3f;1 2f]}];

-1 "\n.htest.ftest:";
.tst.assert_throws["ftest with n=1 group";{.htest.ftest[enlist 1f;1 2 3f]}];

-1 "\n.htest.chisq_gof:";
.tst.assert_throws["chisq_gof length mismatch";{.htest.chisq_gof[10 12 11f;10 10f]}];
.tst.assert_throws["chisq_gof k=1";{.htest.chisq_gof[enlist 10f;enlist 10f]}];

-1 "\n.htest.chisq_ind:";
.tst.assert_throws["chisq_ind too small (1xN)";{.htest.chisq_ind enlist 10 15 20f}];

-1 "\n.htest.anova1:";
.tst.assert_throws["anova1 length mismatch";{.htest.anova1[1 2 3 4f;1 2 3]}];
.tst.assert_throws["anova1 only 1 group";{.htest.anova1[1 2 3f;1 1 1]}];

-1 "\n.htest.cortest:";
.tst.assert_throws["cortest length mismatch";{.htest.cortest[1 2 3f;1 2f]}];
.tst.assert_throws["cortest n<=2";{.htest.cortest[1 2f;3 4f]}];

-1 "\n.htest.proptest:";
.tst.assert_throws["proptest p0<=0";{.htest.proptest[10;100;0f]}];
.tst.assert_throws["proptest p0>=1";{.htest.proptest[10;100;1f]}];
.tst.assert_throws["proptest n=0";{.htest.proptest[10;0;0.5]}];
.tst.assert_throws["proptest x>n";{.htest.proptest[150;100;0.5]}];
.tst.assert_throws["proptest x<0";{.htest.proptest[-1;100;0.5]}];

-1 "\n.htest sanity checks (golden small examples):";

/ One-sample t-test against mean of x: stat ~ 0 and p ~ 1
r1:.htest.ttest1[1 2 3 4 5f;3f];
.tst.assert_approx["ttest1 stat at mean";r1`statistic;0f;1e-12];
.tst.assert_approx["ttest1 p at mean";r1`p_value;1f;1e-12];

/ Two-sample t-test result has the standard 6-key dictionary.
r2:.htest.ttest2[1 2 3 4 5f;2 3 4 5 6f];
.tst.assert["ttest2 returns standard 6-key dict";6=count[r2]];
.tst.assert_equal["ttest2 result keys";key r2;`statistic`df`p_value`method`alternative`ci];

/ Perfect correlation -> p-value ~ 0
r3:.htest.cortest[1 2 3 4 5f;2 4 6 8 10f];
.tst.assert["cortest perfect correlation -> p<1e-10";r3[`p_value]<1e-10];

/ Equal proportions -> z ~ 0, p ~ 1 (p uses pnorm; 5e-7 precision floor)
r4:.htest.proptest[50;100;0.5];
.tst.assert_approx["proptest at p0 stat=0";r4`statistic;0f;1e-12];
.tst.assert_approx["proptest at p0 p=1";r4`p_value;1f;1e-6];
