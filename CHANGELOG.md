# Changelog

All notable changes to qstats are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] — 2026-06-09

First stable release. 106 public functions across 7 namespaces, validated against scipy / statsmodels by 1162 automated tie-out tests.

### Added

#### `.special` — Special functions (Wave 1)
- `lgamma` via Lanczos approximation (g=7, n=9) with reflection for z < 0.5 and recurrence for z < 0.1
- `betainc` via Cephes-style dispatch — power series for small parameters (`b·x ≤ 1 ∧ x ≤ 0.5`), Lentz continued fraction otherwise
- `gammainc` via series expansion for small x, continued fraction for large x

#### `.dist` — 10 probability distributions × {d, p, q, r} (Waves 1, 2, 7)
- Continuous: normal, chi-squared, Student's t, F, uniform, beta, gamma, exponential
- Discrete: binomial, Poisson
- All quantile functions use Halley's method (3rd-order) with Cornish-Fisher / Wilson-Hilferty initial guess, safeguarded with bisection fallback
- Random variates: Box-Muller (normal), Marsaglia-Tsang (gamma), Knuth multiplicative (Poisson), transformation method (others)

#### `.la` — Linear algebra utilities (Wave 1, Wave 6)
- `diag`, `diag_matrix`, `trace`, `det`, `solve`, `crossprod`, `outer`, `eye`, `is_symmetric`
- `cond` — 2-norm condition number via power iteration on `AᵀA` (matches `numpy.linalg.cond` default)
- `power_iter` — dominant-eigenvalue helper
- `eigen_jacobi` — symmetric-matrix eigenvalue solver (used by Durbin-Watson p-value)

#### `.desc` — Descriptive statistics (Wave 3)
- Central tendency: `mean`, `median`, `mode`, `trimmed_mean`
- Dispersion: `var`, `sd`, `iqr`, `mad` (scaled to scipy `scale='normal'`), `cv`
- Shape: `skewness`, `kurtosis` (both type-2, matching R/SAS/SPSS and scipy `bias=False`)
- Order statistics: `quantile` (type-7, R/numpy default), `percentile`
- Correlation: `cor` (Pearson), `spearman`, `kendall` (tau-b)
- Matrices: `cov`, `cov_matrix`, `cor_matrix`
- Tables: `summary`, `freq`, `crosstab`

#### `.htest` — 15 hypothesis tests (Waves 4, 5) with consistent return shape
- Parametric: `ttest1`, `ttest2` (equal-variance pooled), `welch` (Welch-Satterthwaite df), `ttest_paired`, `ftest` (variance ratio), `chisq_gof`, `chisq_ind` (with Yates correction on 2×2), `anova1`, `cortest`, `proptest`
- Nonparametric: `ks` (2-sample, asymptotic Kolmogorov SF), `mannwhitney` (with tie + continuity correction), `wilcoxon` (paired, with tie + continuity correction), `shapiro` (Royston 1992 small-n and large-n branches), `jarque_bera`
- Real confidence intervals on `ttest1`, `ttest2`, `welch`, `ttest_paired`, `ftest`, `cortest`, `proptest`. Others return `(0n; 0n)` where a single-parameter CI is not statistically meaningful.

#### `.lm` + `.diag` — OLS + regression diagnostics (Wave 6)
- `lm.fit[X; y]` — minimal OLS returning a model dictionary
- Influence: `vif`, `cooks_distance`, `leverage`, `dfbetas`, `dffits`
- Heteroskedasticity: `breusch_pagan` (Koenker's studentized form, matches statsmodels default), `white_test`
- Autocorrelation: `durbin_watson` (with p-value via Pan's beta approximation)
- Normality: `jarque_bera` (wrapper over `.htest.jarque_bera` on residuals)
- Model selection: `rsquared` (R² + adjusted), `aic`, `bic` (full statsmodels formula with normal-error log-likelihood)
- Plot support: `residual_data`

### Quality

- **1162 tie-out + unit tests, 0 failures** as of v1.0.0
- Per-function tolerances documented in [README §Numerical Accuracy](README.md#numerical-accuracy)
- Two acknowledged precision floors, both intrinsic to the algorithms used (Shapiro n=4 W; Jacobi solver complexity)
- Test framework's `.tst.assert_approx` rejects null/NaN actuals (q's `0n` sorts below every float, so the original assertion silently passed NaN — fixed in Wave 2)

### Notable algorithm choices (post v0 retrospective)

- **Halley's method for tail quantiles** (Wave 8): replaced naive Newton-Raphson on `qchisq`/`qt` extreme tails. Tightened tolerance from 1e-6 to 1e-10 vs scipy.
- **Safeguarded Newton-Raphson** with bisection fallback used on `qf`, `qbeta`, `qgamma` (Wave 2 / Wave 7): naive NR diverges in tails because the PDF collapses.
- **Yates' continuity correction** on `chisq_ind` 2×2 tables (Wave 4): matches R `chisq.test` and scipy `chi2_contingency` defaults.
- **Pan's beta approximation** for Durbin-Watson (Wave 8): provides p-value via beta-CDF of eigenvalue-derived parameters; tolerance ~1e-2 (intrinsic to Pan).
- **Royston 1992 / Algorithm AS R94** for Shapiro-Wilk (Wave 5 + Wave 8): both n ∈ [12, 5000] log-normal-transform branch and n ∈ [4, 11] γ-transform branch.

### Acknowledged limitations

- `.la.eigen_jacobi` is O(n³); fine to n ≈ 500.
- `.htest.shapiro` n=4 W has a ~2.2e-4 floor vs scipy (single-coefficient denominator collapse — intrinsic to Royston at the smallest sample size).
- Random variates use q's PRNG (not Mersenne Twister), so seed-matched scipy tie-outs aren't possible; correctness is verified via property-based tests at large n.

[1.0.0]: https://github.com/jsogarro/qstats/releases/tag/v1.0.0
