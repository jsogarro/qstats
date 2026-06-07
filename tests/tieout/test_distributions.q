/ test_distributions.q — Tie-out tests for normal distribution
/ Compare against scipy reference values

-1 "\n--- Normal Distribution Tie-Out Tests ---";

/ Load reference data
ref:.j.k raze read0 `:tests/reference/distributions.json;
norm_ref:ref`normal;

/ Test dnorm
-1 "\nTesting .dist.dnorm:";
dnorm_data:norm_ref`dnorm;
dnorm_x:dnorm_data`x;
dnorm_mu:dnorm_data`mu;
dnorm_sigma:dnorm_data`sigma;
dnorm_expected:dnorm_data`y;
dnorm_actual:.dist.dnorm[dnorm_x;dnorm_mu;dnorm_sigma];
{.tst.assert_approx["dnorm(x=",string[x],")";y;z;1e-10]}'[dnorm_x;dnorm_actual;dnorm_expected];

/ Test pnorm (standard normal)
-1 "\nTesting .dist.pnorm (standard):";
pnorm_data:norm_ref`pnorm;
pnorm_x:pnorm_data`x;
pnorm_mu:pnorm_data`mu;
pnorm_sigma:pnorm_data`sigma;
pnorm_expected:pnorm_data`y;
pnorm_actual:.dist.pnorm[pnorm_x;pnorm_mu;pnorm_sigma];
{.tst.assert_approx["pnorm(x=",string[x],")";y;z;5e-7]}'[pnorm_x;pnorm_actual;pnorm_expected];

/ Test pnorm (non-standard)
-1 "\nTesting .dist.pnorm (non-standard):";
pnorm_ns:norm_ref`pnorm_nonstandard;
pnorm_ns_x:pnorm_ns`x;
pnorm_ns_mu:pnorm_ns`mu;
pnorm_ns_sigma:pnorm_ns`sigma;
pnorm_ns_expected:pnorm_ns`y;
pnorm_ns_actual:.dist.pnorm[pnorm_ns_x;pnorm_ns_mu;pnorm_ns_sigma];
{.tst.assert_approx["pnorm(x=",string[x],", mu=5, sigma=2)";y;z;5e-7]}'[pnorm_ns_x;pnorm_ns_actual;pnorm_ns_expected];

/ Test qnorm
-1 "\nTesting .dist.qnorm:";
qnorm_data:norm_ref`qnorm;
qnorm_p:qnorm_data`p;
qnorm_mu:qnorm_data`mu;
qnorm_sigma:qnorm_data`sigma;
qnorm_expected:qnorm_data`y;
qnorm_actual:.dist.qnorm[qnorm_p;qnorm_mu;qnorm_sigma];
{.tst.assert_approx["qnorm(p=",string[x],")";y;z;1e-7]}'[qnorm_p;qnorm_actual;qnorm_expected];

/=============================================================================
/ CHI-SQUARED DISTRIBUTION TESTS
/=============================================================================

-1 "\n--- Chi-Squared Distribution Tie-Out Tests ---";

chisq_ref:ref`chisq;

/ Test dchisq, pchisq, qchisq for each df
test_chisq_df:{[df_key]
  -1 "\nTesting chi-squared with ",string[df_key],":";
  df_data:chisq_ref[df_key];
  df:df_data`df;

  / Test dchisq
  dchisq_data:df_data`dchisq;
  dchisq_x:dchisq_data`x;
  dchisq_expected:dchisq_data`y;
  dchisq_actual:.dist.dchisq[dchisq_x;df];
  {.tst.assert_approx["dchisq(x=",string[x],")";y;z;1e-10]}'[dchisq_x;dchisq_actual;dchisq_expected];

  / Test pchisq
  pchisq_data:df_data`pchisq;
  pchisq_x:pchisq_data`x;
  pchisq_expected:pchisq_data`y;
  pchisq_actual:.dist.pchisq[pchisq_x;df];
  {.tst.assert_approx["pchisq(x=",string[x],")";y;z;1e-10]}'[pchisq_x;pchisq_actual;pchisq_expected];

  / Test qchisq (keep original 1e-6 tolerance for df=1 exact formula)
  qchisq_data:df_data`qchisq;
  qchisq_p:qchisq_data`p;
  qchisq_expected:qchisq_data`y;
  qchisq_actual:.dist.qchisq[qchisq_p;df];
  {.tst.assert_approx["qchisq(p=",string[x],")";y;z;1e-6]}'[qchisq_p;qchisq_actual;qchisq_expected];

  / Test qchisq extreme tails (Wave 8 Halley precision)
  qchisq_ext_data:df_data`qchisq_extreme;
  qchisq_ext_p:qchisq_ext_data`p;
  qchisq_ext_expected:qchisq_ext_data`y;
  qchisq_ext_actual:.dist.qchisq[qchisq_ext_p;df];
  / Use 1e-6 tolerance to accommodate df=1 exact formula precision floor
  {.tst.assert_approx["qchisq_extreme(p=",string[x],")";y;z;1e-6]}'[qchisq_ext_p;qchisq_ext_actual;qchisq_ext_expected];
 };

test_chisq_df each `df_1`df_2`df_5`df_10`df_20;

/=============================================================================
/ STUDENT'S T DISTRIBUTION TESTS
/=============================================================================

-1 "\n--- Student's t Distribution Tie-Out Tests ---";

t_ref:ref`t;

/ Test dt, pt, qt for each df
test_t_df:{[df_key]
  -1 "\nTesting t-distribution with ",string[df_key],":";
  df_data:t_ref[df_key];
  df:df_data`df;

  / Test dt
  dt_data:df_data`dt;
  dt_x:dt_data`x;
  dt_expected:dt_data`y;
  dt_actual:.dist.dt[dt_x;df];
  {.tst.assert_approx["dt(x=",string[x],")";y;z;1e-10]}'[dt_x;dt_actual;dt_expected];

  / Test pt
  pt_data:df_data`pt;
  pt_x:pt_data`x;
  pt_expected:pt_data`y;
  pt_actual:.dist.pt[pt_x;df];
  {.tst.assert_approx["pt(x=",string[x],")";y;z;1e-10]}'[pt_x;pt_actual;pt_expected];

  / Test qt (keep original 1e-6 tolerance)
  qt_data:df_data`qt;
  qt_p:qt_data`p;
  qt_expected:qt_data`y;
  qt_actual:.dist.qt[qt_p;df];
  {.tst.assert_approx["qt(p=",string[x],")";y;z;1e-6]}'[qt_p;qt_actual;qt_expected];

  / Test qt extreme tails (Wave 8 Halley precision)
  qt_ext_data:df_data`qt_extreme;
  qt_ext_p:qt_ext_data`p;
  qt_ext_expected:qt_ext_data`y;
  qt_ext_actual:.dist.qt[qt_ext_p;df];
  / Use 5e-8 tolerance to accommodate betainc precision floor at extreme tails for df=1,2
  {.tst.assert_approx["qt_extreme(p=",string[x],")";y;z;5e-8]}'[qt_ext_p;qt_ext_actual;qt_ext_expected];
 };

test_t_df each `df_1`df_2`df_5`df_10`df_30;

/=============================================================================
/ F DISTRIBUTION TESTS
/=============================================================================

-1 "\n--- F Distribution Tie-Out Tests ---";

f_ref:ref`f;

/ Test df, pf, qf for each df pair
test_f_df:{[df_key]
  -1 "\nTesting F-distribution with ",string[df_key],":";
  df_data:f_ref[df_key];
  df1:df_data`df1;
  df2:df_data`df2;

  / Test df (PDF)
  df_pdf_data:df_data`df;
  df_x:df_pdf_data`x;
  df_expected:df_pdf_data`y;
  df_actual:.dist.df[df_x;df1;df2];
  {.tst.assert_approx["df(x=",string[x],")";y;z;1e-10]}'[df_x;df_actual;df_expected];

  / Test pf
  pf_data:df_data`pf;
  pf_x:pf_data`x;
  pf_expected:pf_data`y;
  pf_actual:.dist.pf[pf_x;df1;df2];
  {.tst.assert_approx["pf(x=",string[x],")";y;z;1e-10]}'[pf_x;pf_actual;pf_expected];

  / Test qf
  qf_data:df_data`qf;
  qf_p:qf_data`p;
  qf_expected:qf_data`y;
  qf_actual:.dist.qf[qf_p;df1;df2];
  / 1e-6 matches qt/qchisq tolerances and the plan's quantile-accuracy spec
  / (1e-6 central, ~4e-6 extreme tails). The residual at p=0.99, df=(1,1)
  / propagates from betainc precision near z=1 (Wave 1 limitation).
  {.tst.assert_approx["qf(p=",string[x],")";y;z;1e-6]}'[qf_p;qf_actual;qf_expected];
 };

test_f_df each `df_1_1`df_5_2`df_5_10`df_10_10`df_20_20;

/=============================================================================
/ UNIFORM DISTRIBUTION TESTS
/=============================================================================

-1 "\n--- Uniform Distribution Tie-Out Tests ---";

unif_ref:ref`uniform;

/ Test dunif (0,1)
-1 "\nTesting .dist.dunif (0,1):";
dunif_data:unif_ref`dunif_01;
dunif_x:dunif_data`x;
dunif_a:dunif_data`a;
dunif_b:dunif_data`b;
dunif_expected:dunif_data`y;
dunif_actual:.dist.dunif[dunif_x;dunif_a;dunif_b];
{.tst.assert_approx["dunif(x=",string[x],")";y;z;1e-10]}'[dunif_x;dunif_actual;dunif_expected];

/ Test punif (0,1)
-1 "\nTesting .dist.punif (0,1):";
punif_data:unif_ref`punif_01;
punif_x:punif_data`x;
punif_a:punif_data`a;
punif_b:punif_data`b;
punif_expected:punif_data`y;
punif_actual:.dist.punif[punif_x;punif_a;punif_b];
{.tst.assert_approx["punif(x=",string[x],")";y;z;1e-10]}'[punif_x;punif_actual;punif_expected];

/ Test qunif (0,1)
-1 "\nTesting .dist.qunif (0,1):";
qunif_data:unif_ref`qunif_01;
qunif_p:qunif_data`p;
qunif_a:qunif_data`a;
qunif_b:qunif_data`b;
qunif_expected:qunif_data`y;
qunif_actual:.dist.qunif[qunif_p;qunif_a;qunif_b];
{.tst.assert_approx["qunif(p=",string[x],")";y;z;1e-10]}'[qunif_p;qunif_actual;qunif_expected];

/ Test dunif (custom range 5,10)
-1 "\nTesting .dist.dunif (5,10):";
dunif_c:unif_ref`dunif_custom;
dunif_c_x:dunif_c`x;
dunif_c_a:dunif_c`a;
dunif_c_b:dunif_c`b;
dunif_c_expected:dunif_c`y;
dunif_c_actual:.dist.dunif[dunif_c_x;dunif_c_a;dunif_c_b];
{.tst.assert_approx["dunif(x=",string[x],", a=5, b=10)";y;z;1e-10]}'[dunif_c_x;dunif_c_actual;dunif_c_expected];
