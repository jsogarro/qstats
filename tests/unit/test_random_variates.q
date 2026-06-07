/ test_random_variates.q -- Statistical property tests for .dist.r* generators
/ Per Wave 2 plan sections 1.4, 2.4, 3.4, 4.4.
/ Approach: count fraction of samples <= theoretical p-th quantile; should ~= p.

-1 "\n--- Random Variate Statistical Tests ---";

N:10000;
TOL:0.015;

.tst.assert_emp_quantile:{[nm;samples;threshold;p]
  / Guard: if the theoretical quantile returned null (q/t/chisq/f Newton-Raphson
  / divergence), report that explicitly rather than masking it as a low emp pct.
  if[null threshold;
    .tst.fail+:1; .tst.failures,:enlist nm;
    -1 "  FAIL: ",nm," (threshold is null — likely q* Newton-Raphson divergence)";
    :()];
  emp:(count samples where samples<=threshold)%count samples;
  .tst.assert_approx[nm;`float$emp;`float$p;TOL]
 };

/ ====================================================================
/ CHI-SQUARED
/ ====================================================================
-1 "\nTesting .dist.rchisq:";

test_rchisq:{[df]
  samples:.dist.rchisq[N;df];
  thresh:.dist.qchisq[0.95;df];
  .tst.assert_emp_quantile["rchisq df=",(string df)," 95th pct";samples;thresh;0.95]
 };
test_rchisq each 1f,2f,5f,10f,20f;

samp_c:.dist.rchisq[N;10f];
.tst.assert_approx["rchisq df=10 mean";avg samp_c;10f;0.2];
.tst.assert_approx["rchisq df=10 var";var samp_c;20f;0.8];

/ ====================================================================
/ STUDENT'S T
/ ====================================================================
-1 "\nTesting .dist.rt:";

test_rt:{[df]
  samples:.dist.rt[N;df];
  thresh:.dist.qt[0.975;df];
  .tst.assert_emp_quantile["rt df=",(string df)," 97.5th pct";samples;thresh;0.975]
 };
test_rt each 1f,2f,5f,10f,30f;

samp_t:.dist.rt[N;5f];
med_t:(asc samp_t)`long$N%2;
.tst.assert_approx["rt df=5 median";med_t;0f;0.1];

/ ====================================================================
/ F
/ ====================================================================
-1 "\nTesting .dist.rf:";

test_rf:{[df1;df2]
  samples:.dist.rf[N;df1;df2];
  thresh:.dist.qf[0.95;df1;df2];
  .tst.assert_emp_quantile["rf df=",(string df1),",",(string df2)," 95th pct";samples;thresh;0.95]
 };
test_rf[1f;1f];
test_rf[5f;2f];
test_rf[5f;10f];
test_rf[10f;10f];
test_rf[20f;20f];

/ ====================================================================
/ UNIFORM
/ ====================================================================
-1 "\nTesting .dist.runif:";

samp_u:.dist.runif[N;0f;1f];
.tst.assert_approx["runif(0,1) mean";avg samp_u;0.5;0.02];
.tst.assert["runif(0,1) min in range";((min samp_u)>=0f) & (min samp_u)<1f];
.tst.assert["runif(0,1) max in range";((max samp_u)>0f) & (max samp_u)<=1f];
.tst.assert_emp_quantile["runif(0,1) 50th pct";samp_u;0.5;0.5];
.tst.assert_emp_quantile["runif(0,1) 25th pct";samp_u;0.25;0.25];

samp_u2:.dist.runif[N;5f;10f];
.tst.assert_approx["runif(5,10) mean";avg samp_u2;7.5;0.1];
.tst.assert["runif(5,10) min in range";((min samp_u2)>=5f) & (min samp_u2)<10f];
.tst.assert["runif(5,10) max in range";((max samp_u2)>5f) & (max samp_u2)<=10f];
.tst.assert_emp_quantile["runif(5,10) 50th pct";samp_u2;7.5;0.5];

samp_u3:.dist.runif[N;-2f;3f];
.tst.assert_approx["runif(-2,3) mean";avg samp_u3;0.5;0.1];
