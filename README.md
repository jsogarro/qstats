# qstats

Statistical testing and diagnostics library for Q/kdb+. Pure q, zero dependencies.

## Why qstats?

The kdb+ ecosystem has no lightweight statistics library. Everyone bridges to Python or writes ad-hoc functions. qstats provides production-quality implementations of the functions you need most — **106 functions across 7 namespaces** covering distributions, descriptive statistics, hypothesis tests, regression diagnostics, and special functions — all in pure q.

Every function is numerically validated against scipy / statsmodels via an automated tie-out suite (**1162 tests, 0 failures**). See [Numerical Accuracy](#numerical-accuracy) for per-function precision targets and known limitations.

## Quick Start

```q
\l src/load.q

/ ---- Distributions: d/p/q/r for normal, t, F, χ², uniform, beta, gamma, binomial, Poisson, exponential ----
.dist.pnorm[1.96; 0; 1]                / 0.975 (CDF)
.dist.qt[0.975; 10]                    / 2.228 (Student's t quantile)
.dist.pchisq[3.84; 1]                  / 0.95  (chi-squared CDF)
.dist.pbeta[enlist 0.5; 2.0; 5.0]      / 0.891 (Beta CDF)
.dist.dpois[3f; 5.0]                   / 0.140 (Poisson PMF)
.dist.pexp[1.0; 2.0]                   / 0.865 (Exponential CDF)

/ ---- Descriptive statistics ----
.desc.mean    1 2 3 4 5f                                  / 3.0
.desc.sd      1 2 3 4 5f                                  / 1.581
.desc.cor    [1 2 3 4 5f; 2 4 6 8 10f]                    / 1.0  (Pearson)
.desc.quantile[1 2 3 4 5 6 7 8 9 10f; 0.75]               / 7.75 (type-7)
.desc.skewness 1 2 3 4 10f                                / 1.515 (type-2)

/ ---- Hypothesis tests ----
.htest.welch    [.dist.rnorm[100;0f;1f]; .dist.rnorm[100;2f;1f]]   / Welch's t-test
.htest.cortest  [1 2 3 4 5f; 2 4 6 8 10f]                          / Pearson correlation test
.htest.ks       [.dist.rnorm[50;0f;1f]; .dist.rnorm[50;1f;1f]]     / 2-sample Kolmogorov-Smirnov
.htest.wilcoxon [before; after]                                    / paired signed-rank

/ ---- OLS + regression diagnostics ----
m: .lm.fit[X; y]                       / fit, returns model dict
.diag.vif             m                / Variance Inflation Factor per predictor
.diag.cooks_distance  m                / Cook's distance per observation
.diag.breusch_pagan   m                / heteroskedasticity LM test
.diag.durbin_watson   m                / autocorrelation test (with p-value)
.diag.aic             m                / Akaike information criterion
```

## Installation

Requires kdb+ 3.5+ (64-bit recommended).

```bash
git clone https://github.com/jsogarro/qstats.git
cd qstats
q src/load.q
```

## Testing

```bash
# Generate scipy / statsmodels reference values (one-time; requires scipy, numpy, statsmodels)
cd tests/reference && python gen_all.py && cd ../..

# Run the full tie-out + unit suite
q tests/run_all.q
# Expected: 1162 passed, 0 failed
```

### Random-Variate Generators

qstats uses a different PRNG than scipy / R, so deterministic tie-out tests are not applicable for `.dist.r*` functions. Correctness is verified via statistical properties: empirical quantile coverage matches theoretical values within ±0.015 at n=10,000 (>99% confidence). This approach aligns with NIST SP 800-22 and provides equivalent guarantees to deterministic tests for continuous distributions.

## Namespaces

| Namespace | Module | Surface |
|-----------|--------|---------|
| `.special` | `src/special.q` | Log-gamma (Lanczos), regularized incomplete beta (Cephes-style dispatch: power series + Lentz CF), regularized incomplete gamma |
| `.dist` | `src/distributions.q` | 10 distributions × {d, p, q, r}: normal, chi-squared, Student's t, F, uniform, beta, gamma, binomial, Poisson, exponential |
| `.la` | `src/linalg.q` | Matrix utilities: diag, trace, det, solve, crossprod, outer, eye, is_symmetric, 2-norm condition number (power iteration), Jacobi eigenvalue solver |
| `.desc` | `src/descriptive.q` | Central tendency, dispersion, shape (skewness/kurtosis type-2), order statistics (quantile type-7), correlation (Pearson/Spearman/Kendall), covariance + cor/cov matrices, summary/freq/crosstab |
| `.htest` | `src/tests.q` | 15 hypothesis tests with consistent return shape: t-tests (1-sample, 2-sample, Welch, paired), F-test (variance), χ² (GoF + independence w/ Yates), one-way ANOVA, Pearson correlation, proportion z-test, KS (2-sample), Mann-Whitney U, Wilcoxon signed-rank, Shapiro-Wilk, Jarque-Bera |
| `.lm` | `src/diagnostics.q` | `lm.fit[X; y]` — minimal OLS, returns model dictionary consumed by `.diag.*` |
| `.diag` | `src/diagnostics.q` | Regression diagnostics: VIF, Cook's distance, leverage, DFBETAS, DFFITS, Durbin-Watson (with p-value via Pan), Breusch-Pagan (Koenker), White's test, R² / adjusted R², AIC, BIC, residual-plot data |

Every `.htest.*` returns the same dictionary shape:

```q
`statistic`df`p_value`method`alternative`ci ! (stat; df; pval; method_str; alt_str; (lo; hi))
```

Six tests (`ttest1`, `ttest2`, `welch`, `ttest_paired`, `cortest`, `proptest`, `ftest`) compute a real CI; the rest return `(0n; 0n)` where a single-parameter CI is not statistically meaningful (omnibus tests, distribution-shape tests).

## Numerical Accuracy

qstats targets parity with scipy / statsmodels reference implementations, validated by an automated tie-out suite (`tests/tieout/`).

| Function class | Tolerance vs reference | Notes |
|---|---|---|
| Special functions (`lgamma`, `betainc`, `gammainc`) | 1e-10 absolute | `betainc` dispatches to power series for small parameters (Cephes gate `b·x ≤ 1 ∧ x ≤ 0.5`) and Lentz CF otherwise |
| PDFs (`d*`) | 1e-10 absolute | Log-domain computation throughout |
| CDFs (`p*`) | 1e-10 absolute | Built on regularized incomplete beta / gamma |
| Quantile functions (`qchisq`, `qt`, `qf`, `qbeta`, `qgamma`) | **1e-10 absolute** including tails | Halley's method (3rd-order) with Cornish-Fisher or Wilson-Hilferty initial guess; safeguarded with bisection fallback |
| Quantile functions (`qnorm`, `qunif`, `qexp`, `qbinom`, `qpois`) | 1e-7 (qnorm) / closed-form exact (others) | |
| Random variates | Property-based (see Random-Variate Generators) | Box-Muller, Marsaglia-Tsang, transformation method |
| Descriptive statistics | 1e-12 absolute | |
| Parametric tests (statistic + p-value + CI) | 1e-10 | `proptest` p-value at 1e-6 (pnorm precision floor) |
| Nonparametric tests | 1e-10 statistic; 1e-6 / 5e-6 p-value (pnorm / KS series truncation) | Shapiro-Wilk W at 1e-4, p-value at 1e-2 (Royston 1992 approximation) |
| Regression diagnostics | 1e-10 / 1e-12 | Durbin-Watson p-value at 1e-2 (Pan's beta approximation) |

### Known precision floors

- **`.htest.shapiro` n=4 W statistic** — fixed ~2.2e-4 floor at the smallest sample size due to single-coefficient denominator collapse in the Royston weight construction; intrinsic to the algorithm, not an implementation defect. n ≥ 5 matches scipy within 7e-5.
- **`.la.eigen_jacobi` performance** — O(n³) Jacobi sweeps. Fine for typical regression problems (n ≤ 500); slower than ideal for very large matrices. Optimization (Lanczos / QR algorithm) is a candidate for a future wave.

## Roadmap

- [x] **Wave 1** — Foundations: special functions, normal distribution, linear algebra
- [x] **Wave 2** — Core distributions: chi-squared, Student's t, F, uniform
- [x] **Wave 3** — Descriptive statistics
- [x] **Wave 4** — Parametric hypothesis tests
- [x] **Wave 5** — Nonparametric tests
- [x] **Wave 6** — OLS + regression diagnostics
- [x] **Wave 7** — Additional distributions: beta, gamma, binomial, Poisson, exponential
- [x] **Wave 8** — Precision floors closed (qchisq/qt Halley, Shapiro small-n, Durbin-Watson p-value) + confidence intervals on all standard tests

### Future work

- Pelican documentation site
- Standard CSV fixtures (mtcars, iris, Boston housing) for richer example-driven tests
- CI/CD via GitHub Actions
- Performance benchmarks (`benchmarks/`)
- Lanczos / QR eigenvalue solver for large-n regression diagnostics

## License

MIT
