/ test_dist_extra_validation.q -- Signal traps for Wave 7 distributions.

-1 "\n--- Additional Distributions Validation ---";

.tst.assert_throws["dbeta with alpha<=0";{.dist.dbeta[0.5;0f;1f]}];
.tst.assert_throws["dbeta with beta<=0";{.dist.dbeta[0.5;1f;-1f]}];
.tst.assert_throws["qbeta with p<=0";{.dist.qbeta_scalar[0f;1f;1f]}];
.tst.assert_throws["qbeta with p>=1";{.dist.qbeta_scalar[1f;1f;1f]}];
.tst.assert_throws["rbeta with n<0";{.dist.rbeta[-1;1f;1f]}];

.tst.assert_throws["dgamma with alpha<=0";{.dist.dgamma[1f;0f;1f]}];
.tst.assert_throws["qgamma with p<=0";{.dist.qgamma_scalar[0f;1f;1f]}];

.tst.assert_throws["dbinom with n<0";{.dist.dbinom[1f;-1;0.5]}];
.tst.assert_throws["dbinom with p>1";{.dist.dbinom[1f;10;1.1]}];

.tst.assert_throws["dpois with lambda<=0";{.dist.dpois[1f;0f]}];

.tst.assert_throws["dexp with rate<=0";{.dist.dexp[1f;0f]}];
.tst.assert_throws["qexp with p>=1";{.dist.qexp[1f;1f]}];

/ Sanity checks on small known values
.tst.assert_approx["dbeta(0.5;1;1) = 1 (uniform)";.dist.dbeta[0.5;1f;1f];1f;1e-12];
.tst.assert_approx["pbeta(0.5;1;1) = 0.5";first .dist.pbeta[enlist 0.5;1f;1f];0.5;1e-12];
.tst.assert_approx["dexp(0;1) = 1";.dist.dexp[0f;1f];1f;1e-12];
.tst.assert_approx["pexp(infinity;1) ~ 1";.dist.pexp[100f;1f];1f;1e-12];
.tst.assert_approx["qexp(0.5;1) = log(2)";first .dist.qexp[enlist 0.5;1f];log 2f;1e-12];
.tst.assert_equal["dpois(0;1) = exp(-1)";.dist.dpois[0f;1f];exp -1f];
.tst.assert_equal["qbinom(1;10;0.5) = 10";first .dist.qbinom[enlist 1f;10;0.5];10];
