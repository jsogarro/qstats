/ test_descriptive_extras.q -- Unit tests for descriptive functions with
/ dict/table outputs (summary, freq, crosstab) plus input-validation checks
/ for the rest. Tieout suite covers numeric correctness; these tests cover
/ structure, edge cases, and signaled errors.

-1 "\n--- Descriptive Statistics: Extras & Validation ---";

/ -----------------------------------------------------------------------------
/ summary -- returns 8-number dictionary
/ -----------------------------------------------------------------------------
-1 "\n.desc.summary:";

s:.desc.summary 1 2 3 4 5f;
.tst.assert_equal["summary keys";key s;`n`mean`sd`min`q25`median`q75`max];
.tst.assert_equal["summary n";s`n;5];
.tst.assert_approx["summary mean";s`mean;3f;1e-12];
.tst.assert_approx["summary sd";s`sd;sqrt 2.5;1e-12];
.tst.assert_approx["summary min";s`min;1f;1e-12];
.tst.assert_approx["summary q25";s`q25;2f;1e-12];
.tst.assert_approx["summary median";s`median;3f;1e-12];
.tst.assert_approx["summary q75";s`q75;4f;1e-12];
.tst.assert_approx["summary max";s`max;5f;1e-12];

/ -----------------------------------------------------------------------------
/ freq -- returns table {val, cnt, prop}
/ -----------------------------------------------------------------------------
-1 "\n.desc.freq:";

f:.desc.freq 1 1 2 2 2 3f;
.tst.assert_equal["freq columns";cols f;`val`cnt`prop];
.tst.assert_equal["freq vals";f`val;1 2 3f];
.tst.assert_equal["freq cnts";f`cnt;2 3 1];
.tst.assert_approx["freq props";f`prop;(2 3 1)%6f;1e-12];
.tst.assert["freq props sum to 1";1f=sum f`prop];

f2:.desc.freq `a`b`a`c`a`b;
.tst.assert_equal["freq symbol vals";f2`val;`a`b`c];
.tst.assert_equal["freq symbol cnts";f2`cnt;3 2 1];

/ -----------------------------------------------------------------------------
/ crosstab -- contingency table as a matrix
/ -----------------------------------------------------------------------------
-1 "\n.desc.crosstab:";

xx:1 1 2 2 3 3;
yy:`a`b`a`b`a`b;
ct:.desc.crosstab[xx;yy];
.tst.assert_equal["crosstab shape (3x2)";(count ct;count ct 0);3 2];
.tst.assert_equal["crosstab values (each cell = 1)";ct;(3 2)#1 1 1 1 1 1];

xx2:1 1 1 2 2 3;
yy2:`a`a`b`a`b`b;
ct2:.desc.crosstab[xx2;yy2];
/ Row 1 (x=1): a=2, b=1. Row 2 (x=2): a=1, b=1. Row 3 (x=3): a=0, b=1.
.tst.assert_equal["crosstab uneven";ct2;(3 2)#2 1 1 1 0 1];

/ -----------------------------------------------------------------------------
/ Validation -- assert invalid arguments signal
/ -----------------------------------------------------------------------------
-1 "\n.desc validation:";

.tst.assert_throws["var with n=1";{.desc.var enlist 3f}];
.tst.assert_throws["sd with n=1";{.desc.sd enlist 3f}];
.tst.assert_throws["skewness with n=2";{.desc.skewness 1 2f}];
.tst.assert_throws["kurtosis with n=3";{.desc.kurtosis 1 2 3f}];
.tst.assert_throws["quantile p<0";{.desc.quantile_scalar[1 2 3f;-0.1]}];
.tst.assert_throws["quantile p>1";{.desc.quantile_scalar[1 2 3f;1.1]}];
.tst.assert_throws["quantile empty";{.desc.quantile_scalar[`float$();0.5]}];
.tst.assert_throws["trimmed_mean trim>=0.5";{.desc.trimmed_mean[1 2 3 4f;0.5]}];
.tst.assert_throws["trimmed_mean trim<0";{.desc.trimmed_mean[1 2 3 4f;-0.1]}];
.tst.assert_throws["cv with mean=0";{.desc.cv -1 0 1f}];
.tst.assert_throws["cov length mismatch";{.desc.cov[1 2 3f;1 2f]}];
.tst.assert_throws["cor length mismatch";{.desc.cor[1 2 3f;1 2f]}];
.tst.assert_throws["spearman length mismatch";{.desc.spearman[1 2 3f;1 2f]}];
.tst.assert_throws["kendall length mismatch";{.desc.kendall[1 2 3f;1 2f]}];
.tst.assert_throws["cov_matrix n=1 row";{.desc.cov_matrix enlist 1 2 3f}];
