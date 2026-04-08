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
