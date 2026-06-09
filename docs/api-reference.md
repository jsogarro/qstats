# qstats API Reference

Complete public-function index for qstats v1.0.0. **106 functions across 7 namespaces.**

Every function is validated against scipy or statsmodels — see [Numerical Accuracy](../README.md#numerical-accuracy) for per-class tolerances.

For runnable examples, see [`docs/examples.q`](examples.q).

---

## `.special` — Special functions

Source: `src/special.q`

| Function | Signature | Description |
|---|---|---|
| `lgamma` | `[z]` | Log-gamma function `ln Γ(z)` via Lanczos approximation (g=7, n=9) with reflection for z < 0.5 and recurrence for z < 0.1 |
| `betainc` | `[x; a; b]` | Regularized incomplete beta `I_x(a, b)` — dispatches to power series for small parameters (`b·x ≤ 1 ∧ x ≤ 0.5`) and Lentz CF otherwise (Cephes-style) |
| `gammainc` | `[a; x]` | Regularized incomplete gamma `P(a, x)` — series expansion for small x, continued fraction for large x |

---

## `.dist` — Probability distributions

Source: `src/distributions.q`. Every distribution exposes four functions following R's naming convention:

- `d<dist>` — PDF (continuous) or PMF (discrete)
- `p<dist>` — CDF
- `q<dist>` — Quantile (inverse CDF)
- `r<dist>` — Random variates

| Distribution | Functions | Parameters | Reference |
|---|---|---|---|
| Normal | `dnorm`, `pnorm`, `qnorm`, `rnorm` | `[x or p or n; mu; sigma]` | scipy.stats.norm |
| Chi-squared | `dchisq`, `pchisq`, `qchisq`, `rchisq` | `[x or p or n; df]` | scipy.stats.chi2 |
| Student's t | `dt`, `pt`, `qt`, `rt` | `[x or p or n; df]` | scipy.stats.t |
| F | `df`, `pf`, `qf`, `rf` | `[x or p or n; df1; df2]` | scipy.stats.f |
| Uniform | `dunif`, `punif`, `qunif`, `runif` | `[x or p or n; a; b]` | scipy.stats.uniform |
| Beta | `dbeta`, `pbeta`, `qbeta`, `rbeta` | `[x or p or n; alpha; beta]` | scipy.stats.beta |
| Gamma | `dgamma`, `pgamma`, `qgamma`, `rgamma` | `[x or p or n; alpha; beta]` (shape, rate) | scipy.stats.gamma |
| Binomial | `dbinom`, `pbinom`, `qbinom`, `rbinom` | `[k or p or num; n; p]` | scipy.stats.binom |
| Poisson | `dpois`, `ppois`, `qpois`, `rpois` | `[k or p or n; lambda]` | scipy.stats.poisson |
| Exponential | `dexp`, `pexp`, `qexp`, `rexp` | `[x or p or n; rate]` | scipy.stats.expon |

**Algorithm notes:**
- Quantile functions `qchisq`, `qt`, `qf`, `qbeta`, `qgamma` use **Halley's method** (3rd-order) with Cornish-Fisher / Wilson-Hilferty initial guess; safeguarded with bisection fallback. Tail precision ≤ 1e-10 vs scipy.
- `qnorm` uses the Beasley-Springer-Moro closed form.
- `qbinom`, `qpois` use binary search on the CDF.
- `rnorm` uses Box-Muller; `rgamma` uses Marsaglia-Tsang (2000); `rbinom` is sum-of-Bernoullis; `rpois` is Knuth's multiplicative algorithm.

---

## `.la` — Linear algebra

Source: `src/linalg.q`

| Function | Signature | Description |
|---|---|---|
| `diag` | `[mat]` | Diagonal of a square matrix |
| `diag_matrix` | `[vec]` | Diagonal matrix from a vector |
| `trace` | `[mat]` | Sum of diagonal entries |
| `det` | `[mat]` | Determinant |
| `solve` | `[A; b]` | Solve `A x = b` via `(inv A) mmu b` |
| `crossprod` | `[X]` | `X' X` |
| `outer` | `[x; y]` | Outer product `x · y'` |
| `eye` | `[n]` | n×n identity |
| `is_symmetric` | `[mat; tol]` | Symmetry check within `tol` |
| `cond` | `[mat]` | 2-norm condition number `σ_max / σ_min` via power iteration on `AᵀA` (matches `numpy.linalg.cond` default) |
| `power_iter` | `[A; maxiter; tol]` | Dominant eigenvalue of `A` via power iteration |
| `eigen_jacobi` | `[A]` | Eigenvalues of a symmetric matrix via Jacobi sweeps |

---

## `.desc` — Descriptive statistics

Source: `src/descriptive.q`

### Central tendency

| Function | Signature | Description |
|---|---|---|
| `mean` | `[x]` | Arithmetic mean |
| `median` | `[x]` | 50th percentile (type-7) |
| `mode` | `[x]` | Most frequent value(s) |
| `trimmed_mean` | `[x; tr]` | Mean after discarding `tr` fraction from each tail |

### Dispersion

| Function | Signature | Description |
|---|---|---|
| `var` | `[x]` | Sample variance (unbiased, n-1 denominator) |
| `sd` | `[x]` | Sample standard deviation |
| `iqr` | `[x]` | Inter-quartile range (Q3 − Q1) |
| `mad` | `[x]` | Median absolute deviation, scaled to match `scipy ... scale='normal'` |
| `cv` | `[x]` | Coefficient of variation (`sd / |mean|`) |

### Shape

| Function | Signature | Description |
|---|---|---|
| `skewness` | `[x]` | Type-2 sample skewness (R/SAS/SPSS, scipy `bias=False`) |
| `kurtosis` | `[x]` | Type-2 excess kurtosis |

### Order statistics

| Function | Signature | Description |
|---|---|---|
| `quantile` | `[x; p]` | Type-7 quantile (R/numpy default); p may be scalar or vector |
| `percentile` | `[x; p]` | Same as `quantile` but with percentage input |

### Correlation + covariance

| Function | Signature | Description |
|---|---|---|
| `cor` | `[x; y]` | Pearson product-moment correlation |
| `spearman` | `[x; y]` | Spearman rank correlation |
| `kendall` | `[x; y]` | Kendall tau-b rank correlation (O(n²)) |
| `cov` | `[x; y]` | Sample covariance (n−1 denominator) |
| `cov_matrix` | `[X]` | Covariance matrix of a data matrix |
| `cor_matrix` | `[X]` | Pearson correlation matrix |

### Summary tables

| Function | Signature | Description |
|---|---|---|
| `summary` | `[x]` | 8-number summary dictionary `{n, mean, sd, min, q25, median, q75, max}` |
| `freq` | `[x]` | Frequency table `(val, cnt, prop)` |
| `crosstab` | `[x; y]` | Contingency matrix of two categorical vectors |

---

## `.htest` — Hypothesis tests

Source: `src/tests.q`. **All tests return the same 6-key dictionary:**

```q
`statistic`df`p_value`method`alternative`ci ! (stat; df; pval; method_str; alt_str; (lo; hi))
```

Six tests compute real confidence intervals; the rest return `(0n; 0n)` where a single-parameter CI is not statistically meaningful.

### Parametric

| Function | Signature | Test | CI |
|---|---|---|---|
| `ttest1` | `[x; mu0]` | One-sample t-test | Yes (t-quantile on mean) |
| `ttest2` | `[x; y]` | Two-sample t-test (equal variance, pooled SD) | Yes (mean difference) |
| `welch` | `[x; y]` | Welch's t-test (unequal variance, Satterthwaite df) | Yes (mean difference) |
| `ttest_paired` | `[x; y]` | Paired t-test (reduces to `ttest1` on differences) | Yes (paired mean diff) |
| `ftest` | `[x; y]` | F-test for variance ratio | Yes (F-quantile-inverted) |
| `chisq_gof` | `[observed; expected]` | Chi-squared goodness-of-fit | No (omnibus) |
| `chisq_ind` | `[table]` | Chi-squared test of independence (Yates' correction on 2×2) | No (omnibus) |
| `anova1` | `[y; groups]` | One-way ANOVA F-test | No (omnibus) |
| `cortest` | `[x; y]` | Pearson correlation t-test | Yes (Fisher z transform) |
| `proptest` | `[x; n; p0]` | One-proportion z-test | Yes (Wilson score CI) |

### Nonparametric

| Function | Signature | Test | CI |
|---|---|---|---|
| `ks` | `[x; y]` | Two-sample Kolmogorov-Smirnov (asymptotic Kolmogorov SF) | No (distribution-shape) |
| `mannwhitney` | `[x; y]` | Mann-Whitney U / Wilcoxon rank-sum (tie + continuity correction) | No (distribution-shape) |
| `wilcoxon` | `[x; y]` | Wilcoxon signed-rank (paired; tie + continuity correction) | No (Hodges-Lehmann out of scope) |
| `shapiro` | `[x]` | Shapiro-Wilk normality (Royston 1992 small-n + large-n branches) | No (omnibus normality) |
| `jarque_bera` | `[x]` | Jarque-Bera normality (biased moments to match scipy) | No (omnibus normality) |

---

## `.lm` — OLS regression

Source: `src/diagnostics.q`

| Function | Signature | Description |
|---|---|---|
| `fit` | `[X; y]` | Minimal OLS: returns a model dictionary `{X, y, beta, residuals, fitted, leverage, n, p, rss, tss, XtX_inv}` consumed by `.diag.*` |

**Note:** `X` must include the intercept column (typically `(n#1f),'X_raw`). The returned model dict is the expected input to every `.diag.*` function.

---

## `.diag` — Regression diagnostics

Source: `src/diagnostics.q`. All accept a model dict returned by `.lm.fit`.

### Influence

| Function | Signature | Returns |
|---|---|---|
| `vif` | `[model]` | Variance Inflation Factor per non-intercept predictor (R `car::vif` convention) |
| `leverage` | `[model]` | Hat-matrix diagonal (accessor; computed during `.lm.fit`) |
| `cooks_distance` | `[model]` | Cook's distance per observation (internally studentized) |
| `dffits` | `[model]` | Standardized change in fitted value when each observation is deleted (externally studentized) |
| `dfbetas` | `[model]` | n × p matrix of standardized coefficient deltas (matches statsmodels `OLSInfluence.dfbetas`) |

### Autocorrelation

| Function | Signature | Returns |
|---|---|---|
| `durbin_watson` | `[model]` | Standard `.htest` dict — DW statistic in [0, 4] plus p-value via Pan's beta approximation (~1e-2 tolerance) |

### Heteroskedasticity

| Function | Signature | Returns |
|---|---|---|
| `breusch_pagan` | `[model]` | Standard `.htest` dict — Koenker's studentized LM = n · R² (matches statsmodels default) |
| `white_test` | `[model]` | Standard `.htest` dict — auxiliary regression includes X, X², and pairwise cross-products of non-constant columns |

### Normality of residuals

| Function | Signature | Returns |
|---|---|---|
| `jarque_bera` | `[model]` | Thin wrapper over `.htest.jarque_bera model``residuals` |

### Model selection

| Function | Signature | Returns |
|---|---|---|
| `rsquared` | `[model]` | Dictionary `{rsq, adj_rsq}` |
| `aic` | `[model]` | Akaike information criterion (full statsmodels formula with normal-error log-likelihood) |
| `bic` | `[model]` | Bayesian information criterion |

### Diagnostic plots

| Function | Signature | Returns |
|---|---|---|
| `residual_data` | `[model]` | Table `(fitted, residuals, std_residuals, leverage, cooks_d)` ready for plotting |

---

## Standard test-result dictionary

Every `.htest.*` and `.diag.*` function that returns a hypothesis-test result uses this 6-key dictionary:

| Key | Type | Notes |
|---|---|---|
| `statistic` | float (or pair for F-based tests) | Test statistic |
| `df` | long, pair, or `0n` | Degrees of freedom (or pair for F/ANOVA, or `0n` when not applicable) |
| `p_value` | float | Two-sided p-value |
| `method` | string | Human-readable test name |
| `alternative` | string | Currently always `"two.sided"` — one-sided variants are planned for v2 |
| `ci` | `(lo; hi)` floats, or `(0n; 0n)` | 95% confidence interval where one exists |

This shape is contractual — programmatic consumers can rely on it.
