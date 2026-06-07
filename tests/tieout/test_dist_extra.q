/ test_dist_extra.q -- Tie-out tests for Wave 7 additional distributions:
/ beta, gamma, binomial, Poisson, exponential. Compares against scipy.stats.

-1 "\n--- Additional Distributions Tie-Out ---";

ref:.j.k raze read0 `:tests/reference/dist_extra.json;

/ ---------- Beta ----------
-1 "\nBeta:";
test_beta:{[k]
  g:ref[`beta;k];
  aa:g`alpha; bb:g`beta;
  -1 "  ",string[k],": a=",(string aa),", b=",string bb;
  d:g`dbeta;
  .tst.assert_approx["dbeta(",string[k],")";`float$.dist.dbeta[d`x;aa;bb];`float$d`y;1e-10];
  pd:g`pbeta;
  .tst.assert_approx["pbeta(",string[k],")";`float$.dist.pbeta[pd`x;aa;bb];`float$pd`y;1e-10];
  qd:g`qbeta;
  .tst.assert_approx["qbeta(",string[k],")";`float$.dist.qbeta[qd`p;aa;bb];`float$qd`y;1e-6]
 };
test_beta each key ref`beta;

/ ---------- Gamma ----------
-1 "\nGamma:";
test_gamma:{[k]
  g:ref[`gamma;k];
  aa:g`alpha; bb:g`beta;
  d:g`dgamma;
  .tst.assert_approx["dgamma(",string[k],")";`float$.dist.dgamma[d`x;aa;bb];`float$d`y;1e-10];
  pd:g`pgamma;
  .tst.assert_approx["pgamma(",string[k],")";`float$.dist.pgamma[pd`x;aa;bb];`float$pd`y;1e-10];
  qd:g`qgamma;
  .tst.assert_approx["qgamma(",string[k],")";`float$.dist.qgamma[qd`p;aa;bb];`float$qd`y;1e-6]
 };
test_gamma each key ref`gamma;

/ ---------- Binomial ----------
-1 "\nBinomial:";
test_binom:{[k]
  g:ref[`binom;k];
  nn:g`n; pp:g`p;
  d:g`dbinom;
  .tst.assert_approx["dbinom(",string[k],")";`float$.dist.dbinom[d`k;nn;pp];`float$d`y;1e-10];
  pd:g`pbinom;
  .tst.assert_approx["pbinom(",string[k],")";`float$.dist.pbinom[pd`k;nn;pp];`float$pd`y;1e-10];
  qd:g`qbinom;
  .tst.assert_equal["qbinom(",string[k],")";`long$.dist.qbinom[qd`prob;nn;pp];`long$qd`y]
 };
test_binom each key ref`binom;

/ ---------- Poisson ----------
-1 "\nPoisson:";
test_pois:{[k]
  g:ref[`poisson;k];
  lam:g`lambda;
  d:g`dpois;
  .tst.assert_approx["dpois(",string[k],")";`float$.dist.dpois[d`k;lam];`float$d`y;1e-10];
  pd:g`ppois;
  .tst.assert_approx["ppois(",string[k],")";`float$.dist.ppois[pd`k;lam];`float$pd`y;1e-10];
  qd:g`qpois;
  .tst.assert_equal["qpois(",string[k],")";`long$.dist.qpois[qd`p;lam];`long$qd`y]
 };
test_pois each key ref`poisson;

/ ---------- Exponential ----------
-1 "\nExponential:";
test_exp:{[k]
  g:ref[`expon;k];
  rate:g`rate;
  d:g`dexp;
  .tst.assert_approx["dexp(",string[k],")";`float$.dist.dexp[d`x;rate];`float$d`y;1e-12];
  pd:g`pexp;
  .tst.assert_approx["pexp(",string[k],")";`float$.dist.pexp[pd`x;rate];`float$pd`y;1e-12];
  qd:g`qexp;
  .tst.assert_approx["qexp(",string[k],")";`float$.dist.qexp[qd`p;rate];`float$qd`y;1e-12]
 };
test_exp each key ref`expon;
