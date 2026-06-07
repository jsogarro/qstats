/ test_diagnostics_validation.q -- Validation + sanity tests for Wave 6.

-1 "\n--- Regression Diagnostics Validation ---";

-1 "\n.lm.fit:";
.tst.assert_throws[".lm.fit X/y length mismatch";{.lm.fit[(2 2)#1 2 3 4f;1 2 3f]}];
.tst.assert_throws[".lm.fit n<=p (singular)";{.lm.fit[(2 2)#1 2 3 4f;1 2f]}];

-1 "\n.diag sanity checks:";

/ Simple regression with known answer: y = 2x + 1 + small deterministic noise.
/ Deterministic noise via sin(x*0.7) avoids relying on \S seed semantics that
/ aren't portable across kdb versions.
N:50;
x_col:`float$1+til N;
y_col:1f+(2f*x_col)+0.1*sin x_col*0.7;
X_design:flip ((N#1f);x_col);
mdl:.lm.fit[X_design;y_col];

.tst.assert["beta length matches p";(count mdl`beta)=mdl`p];
.tst.assert["residuals sum ~ 0 (intercept included)";(abs sum mdl`residuals)<1e-10];
.tst.assert["fitted = X mmu beta";(max abs (X_design mmu mdl`beta)-mdl`fitted)<1e-12];
.tst.assert["leverage sums to p";(abs (mdl`p)-(sum mdl`leverage))<1e-10];

/ rsquared should be very high for tightly fit linear data
r2:.diag.rsquared mdl;
.tst.assert["rsq > 0.999 on perfectly linear data";(r2`rsq)>0.999];
.tst.assert["adj_rsq <= rsq";(r2`adj_rsq)<=r2`rsq];

/ AIC/BIC return floats
.tst.assert["aic is finite";not null .diag.aic mdl];
.tst.assert["bic is finite";not null .diag.bic mdl];
.tst.assert["bic > aic when n > exp(2)";(.diag.bic mdl)>.diag.aic mdl];

/ Durbin-Watson returns a dict with statistic, p_value, etc.
/ DW statistic is in [0, 4]. p_value is in [0, 1].
dw_res:.diag.durbin_watson mdl;
.tst.assert["dw has statistic";`statistic in key dw_res];
.tst.assert["dw has p_value";`p_value in key dw_res];
dw_stat:dw_res`statistic;
dw_pval:dw_res`p_value;
.tst.assert["dw statistic in [0, 4]";((dw_stat)>=0f) and dw_stat<=4f];
.tst.assert["dw p_value in [0, 1]";((dw_pval)>=0f) and dw_pval<=1f];

/ Cook's distance / DFFITS / DFBETAS shapes
.tst.assert["cooks_distance length = n";N=count .diag.cooks_distance mdl];
.tst.assert["dffits length = n";N=count .diag.dffits mdl];
dfb:.diag.dfbetas mdl;
.tst.assert["dfbetas shape = (n, p)";((count dfb)=N) and (count first dfb)=mdl`p];

/ VIF length = p - 1 (excludes intercept)
vifs:.diag.vif mdl;
.tst.assert["vif length = p - 1";(count vifs)=(mdl`p)-1];

/ residual_data row count + columns
rd:.diag.residual_data mdl;
.tst.assert["residual_data has n rows";N=count rd];
.tst.assert_equal["residual_data column order";cols rd;`fitted`residuals`std_residuals`leverage`cooks_d];
