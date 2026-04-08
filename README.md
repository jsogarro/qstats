# qstats

Statistical testing and diagnostics library for Q/kdb+. Pure q, zero dependencies.

## Why qstats?

The kdb+ ecosystem has no lightweight statistics library. Everyone bridges to Python or writes ad-hoc functions. qstats provides production-quality implementations of the functions you need most — distributions, hypothesis tests, descriptive stats, and regression diagnostics — all in pure q.

Every function is numerically validated against scipy/R to at least 6 significant figures.

## Quick Start

```q
\l src/load.q

/ Normal distribution
.dist.pnorm[1.96; 0; 1]           / 0.975 (CDF)
.dist.qnorm[0.975; 0; 1]          / 1.96  (quantile)
.dist.dnorm[0; 0; 1]              / 0.399 (PDF)

/ Special functions
.special.lgamma[5.0]              / 3.178 (log-gamma)
.special.betainc[0.5; 2.0; 3.0]  / regularized incomplete beta
.special.gammainc[2.0; 1.0]      / regularized incomplete gamma
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
| `.dist` | `src/distributions.q` | Normal, t, F, chi-sq distributions |
| `.la` | `src/linalg.q` | Matrix diagonal, trace, det, solve |
| `.desc` | `src/descriptive.q` | Mean, median, sd, skewness, kurtosis, correlation |
| `.htest` | `src/tests.q` | t-tests, chi-squared, KS, ANOVA |
| `.diag` | `src/diagnostics.q` | VIF, Cook's distance, Durbin-Watson |

## Roadmap

- [x] Wave 1: Foundations (special functions, normal distribution, linear algebra)
- [ ] Wave 2: Core distributions (t, F, chi-squared)
- [ ] Wave 3: Descriptive statistics
- [ ] Wave 4: Parametric hypothesis tests
- [ ] Wave 5: Nonparametric tests
- [ ] Wave 6: Regression diagnostics
- [ ] Wave 7: Additional distributions

## License

MIT
