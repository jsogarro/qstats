/ load.q — Main loader for qstats library
/ qstats: Statistical Testing & Diagnostics for Q/kdb+
/ License: MIT
/ https://github.com/jsogarro/qstats

\l src/linalg.q
\l src/special.q
\l src/distributions.q

-1 "";
-1 "qstats loaded. Namespaces: .special .dist .la";
-1 "  .special.lgamma[z]          — log-gamma function";
-1 "  .special.betainc[x;a;b]     — regularized incomplete beta";
-1 "  .special.gammainc[a;x]      — regularized incomplete gamma";
-1 "  .dist.dnorm[x;mu;sigma]     — normal PDF";
-1 "  .dist.pnorm[x;mu;sigma]     — normal CDF";
-1 "  .dist.qnorm[p;mu;sigma]     — normal quantile";
-1 "  .dist.rnorm[n;mu;sigma]     — normal random variates";
-1 "";
