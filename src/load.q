/ load.q — Main loader for qstats library
/ qstats: Statistical Testing & Diagnostics for Q/kdb+
/ License: MIT
/ https://github.com/jsogarro/qstats

\l src/linalg.q
\l src/special.q
\l src/distributions.q

-1 "";
-1 "qstats loaded. Namespaces: .special .dist .la";
-1 "  Special functions:";
-1 "    .special.lgamma[z]          — log-gamma function";
-1 "    .special.betainc[x;a;b]     — regularized incomplete beta";
-1 "    .special.gammainc[a;x]      — regularized incomplete gamma";
-1 "  Distributions (PDF/CDF/quantile/random):";
-1 "    .dist.{d|p|q|r}norm         — normal distribution";
-1 "    .dist.{d|p|q|r}chisq        — chi-squared distribution";
-1 "    .dist.{d|p|q|r}t            — Student's t distribution";
-1 "    .dist.{d|p|q|r}f            — F distribution";
-1 "    .dist.{d|p|q|r}unif         — uniform distribution";
-1 "";
