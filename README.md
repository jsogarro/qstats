# qstats

Statistical testing and diagnostics library for Q/kdb+. Pure q, zero dependencies.

## Why qstats?

The kdb+ ecosystem has no lightweight statistics library. Everyone bridges to Python or writes ad-hoc functions. qstats provides production-quality implementations of the functions you need most — distributions, hypothesis tests, descriptive stats, and regression diagnostics — all in pure q.

Every function is numerically validated against scipy/R. See [Numerical Accuracy](#numerical-accuracy) for per-function precision targets and known limitations.

## Quick Start

```q
\l src/load.q

/ Normal distribution
.dist.pnorm[1.96; 0; 1]           / 0.975 (CDF)
.dist.qnorm[0.975; 0; 1]          / 1.96  (quantile)
.dist.dnorm[0; 0; 1]              / 0.399 (PDF)

/ Chi-squared, t, F, uniform (Wave 2)
.dist.pchisq[3.84; 1]             / 0.95  (CDF)
.dist.qt[0.975; 10]               / 2.228 (Student's t quantile)
.dist.pf[2.5; 5; 10]              / 0.901 (F CDF)
.dist.punif[0.5; 0; 1]            / 0.5   (uniform CDF)

/ Special functions
.special.lgamma[5.0]              / 3.178 (log-gamma)
.special.betainc[0.5; 2.0; 3.0]   / regularized incomplete beta
.special.gammainc[2.0; 1.0]       / regularized incomplete gamma
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
# Generate Python reference values (one-time, requires scipy)
cd tests/reference && python gen_all.py && cd ../..

# Run all tests
q tests/run_all.q
```

## Namespaces

| Namespace | Module | Description |
|-----------|--------|-------------|
| `.special` | `src/special.q` | Log-gamma, incomplete beta/gamma |
| `.dist` | `src/distributions.q` | Normal, chi-sq, t, F, uniform distributions (PDF/CDF/quantile/random) |
| `.la` | `src/linalg.q` | Matrix diagonal, trace, det, solve |
| `.desc` | `src/descriptive.q` | Mean, median, sd, skewness, kurtosis, correlation |
| `.htest` | `src/tests.q` | t-tests, chi-squared, KS, ANOVA |
| `.diag` | `src/diagnostics.q` | VIF, Cook's distance, Durbin-Watson |

## Numerical Accuracy

qstats targets parity with scipy/R reference implementations, validated by an automated tie-out suite (`tests/tieout/`) that compares every function against pre-generated scipy values.

| Function class | Tolerance vs scipy | Notes |
|---|---|---|
| Special functions (`lgamma`, `gammainc`, `betainc`) | 1e-10 absolute | Edge case: `betainc` with `a, b ≤ 0.1` falls back to ~1e-4 — small-parameter series fallback is on the roadmap. |
| PDFs (`dnorm`, `dchisq`, `dt`, `df`, `dunif`) | 1e-10 absolute | Log-domain computation throughout. |
| CDFs (`pnorm`, `pchisq`, `pt`, `pf`, `punif`) | 1e-10 absolute | Built on regularized incomplete beta/gamma. |
| Quantile functions (`qchisq`, `qt`, `qf`) | **1e-6 absolute** in central region | Newton-Raphson inversion. **At extreme tails (p ≤ 0.01 or p ≥ 0.99) accumulated NR error is up to ~4e-6.** scipy/R use higher-order methods (Halley, table-driven asymptotics) and reach ~1e-12 here; matching that is on the roadmap. |
| Quantile functions (`qnorm`, `qunif`) | 1e-7 absolute | Closed-form (Beasley-Springer-Moro / linear). |
| Random variates (`rnorm`, `rchisq`, `rt`, `rf`, `runif`) | Distributional (large-n moments match) | Box-Muller, transformation method. |

### Known precision floors
- **`betainc(0.1, 0.1, x)`** — series convergence near the corner; 9 reference points fail at the 1e-10 tolerance. Tracked for the special-functions revisit.
- **`qchisq` / `qt` at p ∈ {0.01, 0.99}** — NR converges to ~1e-6, not 1e-10. Tolerance widened in the corresponding tieout tests; absolute error is documented per call site.

If you need scipy-grade quantile accuracy at extreme tails, consult scipy directly via `embedPy` — qstats does not yet match it there.

## Roadmap

- [x] Wave 1: Foundations (special functions, normal distribution, linear algebra)
- [x] Wave 2: Core distributions (chi-squared, t, F, uniform)
- [ ] Wave 3: Descriptive statistics
- [ ] Wave 4: Parametric hypothesis tests
- [ ] Wave 5: Nonparametric tests
- [ ] Wave 6: Regression diagnostics
- [ ] Wave 7: Additional distributions

## License

MIT
