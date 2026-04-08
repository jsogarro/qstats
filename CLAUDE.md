# qstats — Statistical Testing & Diagnostics for Q/kdb+

## Project Overview
Pure-q statistics library: distributions, hypothesis tests, descriptive stats, regression diagnostics.
Zero external dependencies. Numerically validated against scipy/R.

## Key Commands
```bash
# Load the library
q src/load.q

# Run all tests (unit + tie-out)
q tests/run_all.q

# Generate reference values from Python (one-time)
cd tests/reference && python gen_all.py

# Build documentation site
cd docs-site && make html && make serve
```

## Namespaces
- `.special` — Special functions (lgamma, betainc, gammainc)
- `.dist` — Probability distributions (dnorm, pnorm, qnorm, rnorm, dt, pt, ...)
- `.desc` — Descriptive statistics (mean, median, sd, skewness, kurtosis, cor, ...)
- `.htest` — Hypothesis tests (ttest1, ttest2, welch, chisq, ks, anova, ...)
- `.diag` — Regression diagnostics (vif, cooks_distance, durbin_watson, ...)
- `.la` — Linear algebra utilities (diag, trace, det, solve, cond, ...)

## Testing Strategy
Every function has two test types:
1. **Unit tests** (`tests/unit/`) — edge cases, types, error handling
2. **Tie-out tests** (`tests/tieout/`) — compare output against Python/scipy reference values

Gate rule: a function is NOT complete until its tie-out test passes.

## Architecture
- `src/special.q` — Foundation: lgamma, betainc, gammainc
- `src/distributions.q` — All distribution functions (depends on special.q)
- `src/descriptive.q` — Descriptive statistics
- `src/tests.q` — Hypothesis tests (depends on distributions.q)
- `src/diagnostics.q` — Regression diagnostics
- `src/linalg.q` — Linear algebra utilities
- `src/load.q` — Loads all modules in dependency order

## q Pitfalls (from avlm-q experience)
- NO blank lines inside multi-line function bodies
- Use `  /` (indented) for comments inside functions, never `/` at column 0
- Avoid reserved words as variables: names, lower, upper, type, string, count, sum, avg, etc.
- Lambda closures don't capture local vars — pass explicitly or use while loops
- Matrix inv returns type 0h — use while loops for diagonal extraction
- Use `?[bool_vec;x;y]` not `$[bool_vec;x;y]` for vector conditionals
