/ load.q — Main loader for qstats library
/ qstats: Statistical Testing & Diagnostics for Q/kdb+
/ License: MIT
/ https://github.com/jsogarro/qstats

\l src/linalg.q
\l src/special.q
\l src/distributions.q
\l src/descriptive.q

-1 "";
-1 "qstats loaded. Namespaces: .special .dist .la .desc";
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
-1 "  Descriptive statistics:";
-1 "    .desc.{mean,median,mode,trimmed_mean}";
-1 "    .desc.{var,sd,iqr,mad,cv}";
-1 "    .desc.{skewness,kurtosis,quantile,percentile}";
-1 "    .desc.{cor,spearman,kendall,cov,cov_matrix,cor_matrix}";
-1 "    .desc.{summary,freq,crosstab}";
-1 "";
