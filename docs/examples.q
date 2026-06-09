/=============================================================================
/ examples.q — qstats cookbook
/ Runnable from the repo root: `q docs/examples.q`
/ Each section is self-contained; copy-paste any block into a q session.
/=============================================================================

\l src/load.q

-1 "";
-1 "================================================================";
-1 " qstats cookbook                                                ";
-1 "================================================================";

/=============================================================================
/ 1. Probability distributions
/=============================================================================

-1 "\n--- 1. Probability distributions ---";

/ Every distribution exposes d (PDF/PMF), p (CDF), q (quantile), r (random).
/ Naming follows R: dnorm / pnorm / qnorm / rnorm, etc.

show "Standard normal 97.5th percentile:                ",string .dist.qnorm[0.975;0f;1f];
show "P(Z <= 1.96) under standard normal:               ",string .dist.pnorm[1.96;0f;1f];
show "Density of N(0,1) at 0:                           ",string .dist.dnorm[0f;0f;1f];

show "Student's t(10) 97.5th percentile (critical):     ",string .dist.qt[0.975;10];
show "Chi-squared(5) 95th percentile (critical):        ",string .dist.qchisq[0.95;5];
show "F(5, 10) 95th percentile (critical):              ",string .dist.qf[0.95;5;10];

/ Beta / gamma / Poisson / binomial / exponential (Wave 7)
show "Beta(2, 5) CDF at 0.3:                            ",string first .dist.pbeta[enlist 0.3;2f;5f];
show "Gamma(shape=2, rate=0.5) density at 1:            ",string first .dist.dgamma[enlist 1f;2f;0.5];
show "Poisson(lambda=5) PMF at k=3:                     ",string .dist.dpois[3f;5f];
show "Binomial(n=10, p=0.5) P(X <= 5):                  ",string first .dist.pbinom[enlist 5f;10;0.5];
show "Exponential(rate=2) CDF at 1:                     ",string .dist.pexp[1f;2f];

/ Random variates — n samples, distribution-specific parameters
sample_norm:.dist.rnorm[100;0f;1f];
show "100 standard-normal draws, sample mean:           ",string avg sample_norm;
show "                          sample sd:              ",string .desc.sd sample_norm;

/=============================================================================
/ 2. Descriptive statistics
/=============================================================================

-1 "\n--- 2. Descriptive statistics ---";

xs:1 2 3 4 5 6 7 8 9 10f;

show "Mean / median / SD:                               ",-3!(.desc.mean xs; .desc.median xs; .desc.sd xs);
show "IQR / MAD / CV:                                   ",-3!(.desc.iqr xs; .desc.mad xs; .desc.cv xs);
show "Skewness (type-2) / Kurtosis (type-2):            ",-3!(.desc.skewness xs; .desc.kurtosis xs);
show "Quantile at p=0.25:                               ",string first .desc.quantile[xs;enlist 0.25];

ys:2 4 6 8 10 12 14 16 18 20f;
show "Pearson / Spearman / Kendall(x, y):               ",-3!(.desc.cor[xs;ys]; .desc.spearman[xs;ys]; .desc.kendall[xs;ys]);

/ Summary dict + frequency table
show "Summary dict (n, mean, sd, quartiles, min, max):";
show .desc.summary xs;

/=============================================================================
/ 3. Parametric hypothesis tests
/=============================================================================

-1 "\n--- 3. Parametric hypothesis tests ---";

/ All `.htest.*` return a 6-key dict: statistic, df, p_value, method, alternative, ci

/ One-sample t-test against mu0
show ".htest.ttest1 (x = 1..5, mu0 = 3 -- centered):";
show .htest.ttest1[1 2 3 4 5f;3f];

/ Two-sample Welch's t-test (preferred over pooled for unequal variance)
g1:.dist.rnorm[50;0f;1f];
g2:.dist.rnorm[50;1f;1f];
show "\n.htest.welch  (N(0,1) vs N(1,1), n=50 each):";
show .htest.welch[g1;g2];

/ Correlation test (Pearson) with Fisher z CI
show "\n.htest.cortest (perfect correlation, n=5):";
show .htest.cortest[1 2 3 4 5f;2 4 6 8 10f];

/ One-proportion z-test with Wilson CI
show "\n.htest.proptest (60 / 100 successes, H0: p = 0.5):";
show .htest.proptest[60;100;0.5];

/=============================================================================
/ 4. Nonparametric hypothesis tests
/=============================================================================

-1 "\n--- 4. Nonparametric hypothesis tests ---";

/ Two-sample Kolmogorov-Smirnov
sa:.dist.rnorm[100;0f;1f];
sb:.dist.rnorm[100;0.5;1f];
show ".htest.ks (N(0,1) vs N(0.5,1), n=100 each):";
show .htest.ks[sa;sb];

/ Mann-Whitney U (Wilcoxon rank-sum)
show "\n.htest.mannwhitney (same samples):";
show .htest.mannwhitney[sa;sb];

/ Wilcoxon signed-rank for paired samples
before:1 2 3 4 5 6 7 8 9 10f;
after:1.5 2.4 3.2 4.8 5.1 6.3 7.5 8.0 9.4 10.7;
show "\n.htest.wilcoxon (paired before/after):";
show .htest.wilcoxon[before;after];

/ Jarque-Bera normality test
show "\n.htest.jarque_bera (100 N(0,1) draws -- should fail to reject):";
show .htest.jarque_bera .dist.rnorm[100;0f;1f];

/ Shapiro-Wilk normality test
show "\n.htest.shapiro (30 N(0,1) draws -- should fail to reject):";
show .htest.shapiro .dist.rnorm[30;0f;1f];

/=============================================================================
/ 5. OLS regression + diagnostics
/=============================================================================

-1 "\n--- 5. OLS regression + diagnostics ---";

/ Build a small regression dataset: y = 2 + 1.5*x1 - 0.7*x2 + deterministic noise.
/ x1 and x2 are non-collinear (different shapes) so the design matrix is well-conditioned.
nn:50;
x1:`float$1+til nn;
x2:0.5*sin x1*0.3;
y_true:2f+(1.5*x1)-0.7*x2;
y_noisy:y_true+0.5*sin x1*0.7;

/ Design matrix includes intercept (column of 1s)
X_des:flip ((nn#1f); x1; x2);

/ Fit
mdl:.lm.fit[X_des; y_noisy];
show "Fitted coefficients (intercept, x1, x2):          ",-3!mdl`beta;

/ Diagnostics
r2:.diag.rsquared mdl;
show "R^2 / adjusted R^2:                              ",-3!(r2`rsq; r2`adj_rsq);

show "Variance Inflation Factor (per predictor):       ",-3!.diag.vif mdl;
show "AIC / BIC:                                        ",-3!(.diag.aic mdl; .diag.bic mdl);

show "\nDurbin-Watson test:";
show .diag.durbin_watson mdl;

show "\nBreusch-Pagan (heteroskedasticity) test:";
show .diag.breusch_pagan mdl;

show "\nJarque-Bera on residuals (normality):";
show .diag.jarque_bera mdl;

show "\nDiagnostic plot data (first 5 rows):";
show 5#.diag.residual_data mdl;

/=============================================================================
/ 6. Linear algebra utilities
/=============================================================================

-1 "\n--- 6. Linear algebra utilities ---";

A:(3 3)#1 2 3 4 5 6 7 8 10f;
show "A = ";
show A;
show "diag(A):                                          ",-3!.la.diag A;
show "trace(A):                                         ",string .la.trace A;
show "det(A):                                           ",string .la.det A;
show "cond(A) (2-norm via power iteration):             ",string .la.cond A;

b:1 2 3f;
show "Solve A x = b:                                    ",-3!.la.solve[A;b];

/=============================================================================
/ 7. Special functions (low-level numerics)
/=============================================================================

-1 "\n--- 7. Special functions ---";

show "lgamma(5):                                        ",string .special.lgamma 5f;
show "lgamma(0.5) = ln(sqrt(pi)):                       ",string .special.lgamma 0.5;
show "betainc(0.3, 2, 5):                               ",string .special.betainc[0.3;2f;5f];
show "gammainc(3, 2):                                   ",string .special.gammainc[3f;2f];

-1 "\n================================================================";
-1 " End of cookbook                                                ";
-1 "================================================================";

exit 0;
