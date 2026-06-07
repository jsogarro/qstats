/=============================================================================
/ diagnostics.q -- OLS regression diagnostics
/ Namespaces: .lm (minimal OLS fit), .diag (diagnostics)
/ Wave 6 of qstats master plan.
/ References: statsmodels (matches statsmodels.OLS().fit() and the
/ OLSInfluence diagnostics).
/=============================================================================

.diag.validate:{[check;msg]
  if[not check;'`$"invalid_arg: ",msg];
 };

/=============================================================================
/ Minimal OLS infrastructure (.lm namespace)
/=============================================================================

/ Fit OLS regression y ~ X. X must include the intercept column (if needed).
/ Returns a model dictionary with:
/   X, y, beta, residuals, fitted, leverage, n, p, rss, tss
.lm.fit:{[X;y]
  n:count y;
  p:count first X;
  .diag.validate[n=count X;"X and y must have the same number of rows"];
  .diag.validate[n>p;"more observations than parameters required"];
  XtX:.la.crossprod X;
  Xty:(flip X) mmu y;
  / beta = (X'X)^-1 X'y
  XtX_inv:"f"$inv XtX;
  beta:XtX_inv mmu Xty;
  fitted:X mmu beta;
  resid:y-fitted;
  / Hat matrix diagonal: h_i = x_i^T (X^T X)^-1 x_i.
  / Compute via diagonal of X @ XtX_inv @ X^T without materialising the full
  / n x n hat matrix.
  H_diag:{[XtX_inv;xi] sum xi*XtX_inv mmu xi}[XtX_inv;] each X;
  rss:sum resid*resid;
  tss:sum (y-avg y) xexp 2;
  `X`y`beta`residuals`fitted`leverage`n`p`rss`tss`XtX_inv!(X;y;beta;resid;fitted;H_diag;n;p;rss;tss;XtX_inv)
 };

/=============================================================================
/ 1. Multicollinearity Diagnostics
/=============================================================================

/ Variance Inflation Factor for each predictor (excluding intercept).
/ Assumes the first column of X is the intercept (all 1s).
/ VIF_j = 1 / (1 - R_j^2) where R_j^2 comes from regressing X[;j] on the
/ other predictors (and the intercept).
.diag.vif:{[model]
  X:model`X;
  p:model`p;
  / For each predictor j in 1..p-1 (skip intercept at index 0):
  vif_one:{[X;j]
    n:count X;
    pcols:count first X;
    Xj:X[;j];
    other_cols:where (til pcols)<>j;
    Xo:X[;other_cols];
    / Regress Xj ~ Xo
    XoTXo_inv:"f"$inv .la.crossprod Xo;
    b:XoTXo_inv mmu (flip Xo) mmu Xj;
    fit:Xo mmu b;
    rss_j:sum (Xj-fit) xexp 2;
    tss_j:sum (Xj-avg Xj) xexp 2;
    rsq:1f-rss_j%tss_j;
    1f%1f-rsq
   };
  vif_one[X;] each 1+til p-1
 };

/=============================================================================
/ 2. Influential Observations
/=============================================================================

/ Cook's distance per observation. Uses internal (in-sample) MSE:
/   D_i = (r_i^2 / p) * h_i / (1 - h_i)
/ where r_i is the internally studentized residual:
/   r_i = e_i / sqrt(MSE * (1 - h_i))
.diag.cooks_distance:{[model]
  resid:model`residuals;
  h:model`leverage;
  p:model`p;
  n:model`n;
  mse:(model`rss)%n-p;
  std_resid:resid%sqrt mse*1f-h;
  ((std_resid*std_resid)%p)*h%1f-h
 };

/ Leverage (diagonal of the hat matrix). The model already carries it; this
/ accessor exists so callers don't have to know the dict key.
.diag.leverage:{[model] model`leverage};

/ DFFITS per observation: r_studentized_external * sqrt(h_i / (1 - h_i)).
/ Uses externally studentized residuals (MSE without observation i):
/   sigma_(i)^2 = ((n-p)*MSE - e_i^2/(1-h_i)) / (n-p-1)
.diag.dffits:{[model]
  resid:model`residuals;
  h:model`leverage;
  p:model`p;
  n:model`n;
  rss:model`rss;
  mse:rss%n-p;
  / External studentized variance per observation
  ei2_over_1mh:(resid*resid)%1f-h;
  sigma_ext_sq:((n-p)*mse)-ei2_over_1mh;
  sigma_ext_sq%:(n-p)-1f;
  r_ext:resid%sqrt sigma_ext_sq*1f-h;
  r_ext*sqrt h%1f-h
 };

/ DFBETAS: n x p matrix. For observation i and coefficient j:
/   DFBETAS_{i,j} = ((X^T X)^-1 x_i)_j * e_i / ((1 - h_i) * sigma_(i) * sqrt((X^T X)^-1_{jj}))
/ Matches statsmodels OLSInfluence.dfbetas (externally studentized).
.diag.dfbetas:{[model]
  X:model`X;
  resid:model`residuals;
  h:model`leverage;
  p:model`p;
  n:model`n;
  XtX_inv:model`XtX_inv;
  mse:(model`rss)%n-p;
  ei2_over_1mh:(resid*resid)%1f-h;
  sigma_ext_sq:(((n-p)*mse)-ei2_over_1mh)%(n-p)-1f;
  sigma_ext:sqrt sigma_ext_sq;
  XtX_inv_diag_sqrt:sqrt .la.diag XtX_inv;
  / Per-observation: delta_beta_i = (XtX_inv mmu x_i) * (e_i / (1 - h_i))
  / Then DFBETAS_{i,j} = delta_beta_i[j] / (sigma_ext[i] * sqrt(XtX_inv[j,j]))
  one_row:{[XtX_inv;XtX_inv_diag_sqrt;X;resid;h;sigma_ext;i]
    x_i:X[i];
    e_i:resid[i];
    h_i:h[i];
    sig_i:sigma_ext[i];
    delta_b:(XtX_inv mmu x_i)*e_i%1f-h_i;
    delta_b%sig_i*XtX_inv_diag_sqrt
   };
  one_row[XtX_inv;XtX_inv_diag_sqrt;X;resid;h;sigma_ext;] each til n
 };

/=============================================================================
/ 3. Autocorrelation
/=============================================================================

/ Durbin-Watson statistic for first-order autocorrelation in residuals.
/ DW in [0, 4]; DW=2 means no autocorrelation. p-value is not computed
/ (requires Durbin-Watson tables / Imhof distribution).
.diag.durbin_watson:{[model]
  e:model`residuals;
  numer:sum (1_e)-(-1_e) xexp 2;
  / Right-assoc trap: (1_e)-(-1_e) parses fine since both are parenthesized
  / lists. xexp 2 squares the resulting differences.
  numer:sum ((1_e)-(-1_e)) xexp 2;
  denom:sum e*e;
  numer%denom
 };

/=============================================================================
/ 4. Heteroskedasticity Tests
/=============================================================================

/ Breusch-Pagan LM test, Koenker's studentized form.
/ Regress e^2 on X; LM = n * R^2 ~ chi^2(p-1).
/ Matches statsmodels.stats.diagnostic.het_breuschpagan defaults when X
/ already includes the intercept column. The earlier impl used ESS/2 (the
/ original BP-under-normality formula), which agrees with Koenker only up
/ to a constant of proportionality -- they differ in general.
.diag.breusch_pagan:{[model]
  X:model`X;
  e:model`residuals;
  n:model`n;
  p:model`p;
  u:e*e;
  XtX_inv:"f"$inv .la.crossprod X;
  b_aux:XtX_inv mmu (flip X) mmu u;
  fit_aux:X mmu b_aux;
  rss_aux:sum (u-fit_aux) xexp 2;
  tss_aux:sum (u-avg u) xexp 2;
  rsq:1f-rss_aux%tss_aux;
  lm:n*rsq;
  df:p-1;
  pval:1f-.dist.pchisq[lm;df];
  `statistic`df`p_value`method`alternative`ci!(lm;df;pval;"Breusch-Pagan (Koenker) test";"two.sided";(0n;0n))
 };

/ White's test: auxiliary regression includes original X plus squares and
/ cross-products. LM = n * R^2 ~ chi^2(df) where df = (auxiliary cols - 1).
/ Constructs X_aug = [X, X*X, X cross products of non-constant columns].
.diag.white_test:{[model]
  X:model`X;
  e:model`residuals;
  n:model`n;
  p:model`p;
  / Identify non-intercept columns (those with non-zero variance)
  colidx:til p;
  is_const:{(min x)=max x} each flip X;
  nonconst:colidx where not is_const;
  / Squared columns (only for non-constant cols)
  sq_cols:flip {x*x} each flip X[;nonconst];
  / Pairwise cross products of non-constant columns
  pairs:();
  cn:count nonconst;
  if[cn>1;
    i:0;
    while[i<cn-1;
      j:i+1;
      while[j<cn;
        pairs,:enlist (X[;nonconst[i]])*X[;nonconst[j]];
        j+:1];
      i+:1];
   ];
  X_aug:$[count pairs;
    X,'sq_cols,'flip pairs;
    X,'sq_cols];
  u2:e*e;
  XtX_aug_inv:"f"$inv .la.crossprod X_aug;
  b_aux:XtX_aug_inv mmu (flip X_aug) mmu u2;
  fit_aux:X_aug mmu b_aux;
  rss_aux:sum (u2-fit_aux) xexp 2;
  tss_aux:sum (u2-avg u2) xexp 2;
  rsq:1f-rss_aux%tss_aux;
  lm:n*rsq;
  df:(count first X_aug)-1;
  pval:1f-.dist.pchisq[lm;df];
  `statistic`df`p_value`method`alternative`ci!(lm;df;pval;"White's test for heteroskedasticity";"two.sided";(0n;0n))
 };

/=============================================================================
/ 5. Normality of Residuals
/=============================================================================

/ Jarque-Bera on residuals -- thin wrapper over Wave 5's .htest.jarque_bera.
.diag.jarque_bera:{[model] .htest.jarque_bera model`residuals};

/=============================================================================
/ 6. Model Selection Criteria
/=============================================================================

/ R^2 and adjusted R^2.
.diag.rsquared:{[model]
  rss:model`rss;
  tss:model`tss;
  n:model`n;
  p:model`p;
  r2:1f-rss%tss;
  adj:1f-((1f-r2)*n-1f)%n-p;
  `rsq`adj_rsq!(r2;adj)
 };

/ Log-likelihood for OLS under normal errors. Used by AIC/BIC.
/   ll = -n/2 * (log(2*pi) + 1 + log(RSS / n))
.diag._loglike:{[model]
  n:model`n;
  rss:model`rss;
  / log(2*pi) ~= 1.8378770664093453
  c2pi:log 2f*acos -1f;
  (neg n%2f)*(c2pi+1f)+log rss%n
 };

/ AIC matching statsmodels: AIC = -2 * ll + 2 * p, where p is the parameter
/ count (regressors including intercept).
.diag.aic:{[model]
  ll:.diag._loglike model;
  (neg 2f*ll)+2f*model`p
 };

/ BIC matching statsmodels: BIC = -2 * ll + p * log(n).
.diag.bic:{[model]
  ll:.diag._loglike model;
  (neg 2f*ll)+(model`p)*log model`n
 };

/=============================================================================
/ 7. Residual Plot Data
/=============================================================================

/ Convenience table with fitted, residuals, internally studentized residuals,
/ leverage, and Cook's distance.
.diag.residual_data:{[model]
  fitted:model`fitted;
  resid:model`residuals;
  h:model`leverage;
  n:model`n;
  p:model`p;
  mse:(model`rss)%n-p;
  std_resid:resid%sqrt mse*1f-h;
  cooks:.diag.cooks_distance model;
  ([] fitted:fitted; residuals:resid; std_residuals:std_resid; leverage:h; cooks_d:cooks)
 };
