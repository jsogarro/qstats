/ test_diagnostics.q -- Tie-out tests for Wave 6 regression diagnostics.
/ Compares .lm.fit + .diag.* against statsmodels reference values in
/ tests/reference/diagnostics.json.
/ Tolerances are 1e-10 unless noted; AIC/BIC use exact statsmodels formulas.

-1 "\n--- Regression Diagnostics Tie-Out ---";

ref:.j.k raze read0 `:tests/reference/diagnostics.json;
ds:ref`dataset;

/ Cast nested lists to flat float matrices (JSON gives general lists).
X:"f"$ds`X;
y:"f"$ds`y;

model:.lm.fit[X;y];

/ ----- OLS fit itself -----
-1 "\n.lm.fit:";
.tst.assert_approx["beta";`float$model`beta;`float$ds`beta;1e-10];
.tst.assert_approx["residuals";`float$model`residuals;`float$ds`residuals;1e-10];
.tst.assert_approx["fitted";`float$model`fitted;`float$ds`fitted;1e-10];
.tst.assert_approx["leverage";`float$model`leverage;`float$ds`leverage;1e-10];
.tst.assert_approx["rss";`float$model`rss;`float$ds`rss;1e-8];
.tst.assert_approx["tss";`float$model`tss;`float$ds`tss;1e-8];

/ ----- VIF -----
-1 "\n.diag.vif:";
.tst.assert_approx["vif";`float$.diag.vif model;`float$ref`vif;1e-10];

/ ----- Cook's distance, DFFITS, DFBETAS -----
-1 "\n.diag.{cooks_distance, dffits, dfbetas}:";
.tst.assert_approx["cooks_distance";`float$.diag.cooks_distance model;`float$ref`cooks_distance;1e-10];
.tst.assert_approx["dffits";`float$.diag.dffits model;`float$ref`dffits;1e-10];
dfb_q:`float$.diag.dfbetas model;
dfb_ref:`float$ref`dfbetas;
.tst.assert_approx["dfbetas";`float$raze dfb_q;`float$raze dfb_ref;1e-10];

/ ----- Durbin-Watson -----
-1 "\n.diag.durbin_watson:";
.tst.assert_approx["durbin_watson";`float$.diag.durbin_watson model;`float$ref`durbin_watson;1e-10];

/ ----- Breusch-Pagan -----
-1 "\n.diag.breusch_pagan:";
bp_res:.diag.breusch_pagan model;
bp_ref:ref`breusch_pagan;
.tst.assert_approx["breusch_pagan statistic";`float$bp_res`statistic;`float$bp_ref`statistic;1e-10];
.tst.assert_approx["breusch_pagan p_value";`float$bp_res`p_value;`float$bp_ref`p_value;1e-10];
.tst.assert_approx["breusch_pagan df";`float$bp_res`df;`float$bp_ref`df;1e-10];

/ ----- White's test -----
-1 "\n.diag.white_test:";
wh_res:.diag.white_test model;
wh_ref:ref`white_test;
.tst.assert_approx["white_test statistic";`float$wh_res`statistic;`float$wh_ref`statistic;1e-10];
.tst.assert_approx["white_test p_value";`float$wh_res`p_value;`float$wh_ref`p_value;1e-10];

/ ----- R^2 / adjusted R^2 -----
-1 "\n.diag.rsquared:";
r2_res:.diag.rsquared model;
r2_ref:ref`rsquared;
.tst.assert_approx["rsq";`float$r2_res`rsq;`float$r2_ref`rsq;1e-12];
.tst.assert_approx["adj_rsq";`float$r2_res`adj_rsq;`float$r2_ref`adj_rsq;1e-12];

/ ----- AIC / BIC -----
-1 "\n.diag.aic, .diag.bic:";
.tst.assert_approx["aic";`float$.diag.aic model;`float$ref`aic;1e-8];
.tst.assert_approx["bic";`float$.diag.bic model;`float$ref`bic;1e-8];

/ ----- Jarque-Bera on residuals -----
-1 "\n.diag.jarque_bera:";
jb_res:.diag.jarque_bera model;
jb_ref:ref`jarque_bera;
.tst.assert_approx["jarque_bera statistic";`float$jb_res`statistic;`float$jb_ref`statistic;1e-10];
.tst.assert_approx["jarque_bera p_value";`float$jb_res`p_value;`float$jb_ref`p_value;1e-10];

/ ----- residual_data (structure check) -----
-1 "\n.diag.residual_data:";
rd:.diag.residual_data model;
.tst.assert_equal["residual_data columns";cols rd;`fitted`residuals`std_residuals`leverage`cooks_d];
.tst.assert["residual_data row count";(ds`n)=count rd];
