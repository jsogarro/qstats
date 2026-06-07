#!/usr/bin/env python3
"""Generate reference values for qstats tie-out tests.

Uses scipy/numpy to compute ground truth values that q implementations
must match. Run from the tests/reference/ directory:

    cd tests/reference && python gen_all.py

This regenerates ALL reference JSON files.
"""
import json
import numpy as np
import scipy.special as sp
import scipy.stats as st


def gen_special():
    """Generate reference values for special functions."""
    data = {}

    # Log-gamma. Cover (0.1, 0.5) explicitly: the prior reference set jumped
    # 0.01 -> 0.5, leaving the reflection path (used for 0.1 < z < 0.5)
    # untested. That gap masked the right-associative subtraction bug in the
    # reflection formula until the Wave 3 betainc(a=b=0.1) failures forced a
    # closer look.
    x_lgamma = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 10.0, 50.0, 100.0,
                0.01, 0.001, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49]
    data["lgamma"] = {
        "inputs": {"x": x_lgamma},
        "outputs": {"y": [float(sp.gammaln(x)) for x in x_lgamma]}
    }

    # Regularized incomplete beta I_x(a, b)
    x_beta = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    ab_pairs = [(1.0, 1.0), (2.0, 5.0), (5.0, 2.0), (0.5, 0.5), (10.0, 10.0), (0.1, 0.1)]
    beta_tests = []
    for a, b in ab_pairs:
        vals = [float(sp.betainc(a, b, x)) for x in x_beta]
        beta_tests.append({"a": a, "b": b, "x": x_beta, "y": vals})
    data["betainc"] = beta_tests

    # Regularized incomplete gamma P(a, x)
    x_gamma = [0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0]
    a_vals = [0.5, 1.0, 2.0, 3.0, 5.0, 10.0]
    gamma_tests = []
    for a in a_vals:
        vals = [float(sp.gammainc(a, x)) for x in x_gamma]
        gamma_tests.append({"a": a, "x": x_gamma, "y": vals})
    data["gammainc"] = gamma_tests

    with open("special.json", "w") as f:
        json.dump(data, f, indent=2)
    print(f"  special.json: {len(x_lgamma)} lgamma + {len(beta_tests)} betainc + {len(gamma_tests)} gammainc")


def gen_distributions():
    """Generate reference values for distribution functions."""
    data = {}

    # Normal distribution
    x_norm = [-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
    p_norm = [0.001, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.999]

    data["normal"] = {
        "dnorm": {
            "x": x_norm,
            "mu": 0.0, "sigma": 1.0,
            "y": [float(st.norm.pdf(x, 0, 1)) for x in x_norm]
        },
        "pnorm": {
            "x": x_norm,
            "mu": 0.0, "sigma": 1.0,
            "y": [float(st.norm.cdf(x, 0, 1)) for x in x_norm]
        },
        "qnorm": {
            "p": p_norm,
            "mu": 0.0, "sigma": 1.0,
            "y": [float(st.norm.ppf(p, 0, 1)) for p in p_norm]
        },
        # Non-standard normal
        "pnorm_nonstandard": {
            "x": [0.0, 5.0, 7.0, 10.0],
            "mu": 5.0, "sigma": 2.0,
            "y": [float(st.norm.cdf(x, 5, 2)) for x in [0.0, 5.0, 7.0, 10.0]]
        }
    }

    # Chi-squared distribution
    x_chisq = [0.0, 0.5, 1.0, 2.0, 3.84, 5.0, 10.0, 20.0]
    p_chisq = [0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99]
    df_chisq = [1.0, 2.0, 5.0, 10.0, 20.0]
    # Extreme tail probabilities for Wave 8 Halley precision test
    p_extreme = [0.0001, 0.001, 0.005, 0.01, 0.05, 0.95, 0.99, 0.995, 0.999, 0.9999]

    data["chisq"] = {}
    for df in df_chisq:
        data["chisq"][f"df_{int(df)}"] = {
            "df": df,
            "dchisq": {"x": x_chisq, "y": [float(st.chi2.pdf(x, df)) for x in x_chisq]},
            "pchisq": {"x": x_chisq, "y": [float(st.chi2.cdf(x, df)) for x in x_chisq]},
            "qchisq": {"p": p_chisq, "y": [float(st.chi2.ppf(p, df)) for p in p_chisq]},
            "qchisq_extreme": {"p": p_extreme, "y": [float(st.chi2.ppf(p, df)) for p in p_extreme]}
        }

    # Student's t distribution
    x_t = [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
    p_t = [0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99]
    df_t = [1.0, 2.0, 5.0, 10.0, 30.0]

    data["t"] = {}
    for df in df_t:
        data["t"][f"df_{int(df)}"] = {
            "df": df,
            "dt": {"x": x_t, "y": [float(st.t.pdf(x, df)) for x in x_t]},
            "pt": {"x": x_t, "y": [float(st.t.cdf(x, df)) for x in x_t]},
            "qt": {"p": p_t, "y": [float(st.t.ppf(p, df)) for p in p_t]},
            "qt_extreme": {"p": p_extreme, "y": [float(st.t.ppf(p, df)) for p in p_extreme]}
        }

    # F distribution
    x_f = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 10.0]
    p_f = [0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99]
    df_pairs = [(1, 1), (5, 2), (5, 10), (10, 10), (20, 20)]

    data["f"] = {}
    for df1, df2 in df_pairs:
        data["f"][f"df_{df1}_{df2}"] = {
            "df1": df1, "df2": df2,
            "df": {"x": x_f, "y": [float(st.f.pdf(x, df1, df2)) for x in x_f]},
            "pf": {"x": x_f, "y": [float(st.f.cdf(x, df1, df2)) for x in x_f]},
            "qf": {"p": p_f, "y": [float(st.f.ppf(p, df1, df2)) for p in p_f]}
        }

    # Uniform distribution
    x_unif = [-1.0, 0.0, 0.25, 0.5, 0.75, 1.0, 2.0]
    p_unif = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]

    data["uniform"] = {
        "dunif_01": {
            "x": x_unif, "a": 0.0, "b": 1.0,
            "y": [float(st.uniform.pdf(x, 0, 1)) for x in x_unif]
        },
        "punif_01": {
            "x": x_unif, "a": 0.0, "b": 1.0,
            "y": [float(st.uniform.cdf(x, 0, 1)) for x in x_unif]
        },
        "qunif_01": {
            "p": p_unif, "a": 0.0, "b": 1.0,
            "y": [float(st.uniform.ppf(p, 0, 1)) for p in p_unif]
        },
        "dunif_custom": {
            "x": [0.0, 5.0, 7.5, 10.0, 15.0], "a": 5.0, "b": 10.0,
            "y": [float(st.uniform.pdf(x, 5, 5)) for x in [0.0, 5.0, 7.5, 10.0, 15.0]]  # scale=5 means b-a
        }
    }

    with open("distributions.json", "w") as f:
        json.dump(data, f, indent=2)
    print(f"  distributions.json: normal, chisq, t, F, uniform")


def gen_linalg():
    """Generate reference values for linear algebra functions."""
    np.random.seed(42)
    A = np.random.randn(4, 4)
    A = A @ A.T  # Make symmetric positive definite

    data = {
        "matrix_4x4": {
            "A": A.tolist(),
            "diag": np.diag(A).tolist(),
            "trace": float(np.trace(A)),
            "det": float(np.linalg.det(A)),
            "inv": np.linalg.inv(A).tolist(),
            "cond": float(np.linalg.cond(A))
        }
    }

    # Simple 3x3 for solve
    B = np.array([[2.0, 1.0, -1.0], [-3.0, -1.0, 2.0], [-2.0, 1.0, 2.0]])
    b = np.array([8.0, -11.0, -3.0])
    data["solve_3x3"] = {
        "A": B.tolist(),
        "b": b.tolist(),
        "x": np.linalg.solve(B, b).tolist()
    }

    with open("linalg.json", "w") as f:
        json.dump(data, f, indent=2)
    print(f"  linalg.json: 4x4 matrix ops + 3x3 solve")


def gen_descriptive():
    """Generate reference values for descriptive statistics (Wave 3)."""
    rng = np.random.default_rng(42)
    x_small = [1.0, 2.0, 3.0, 4.0, 5.0]
    x_med = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    x_mixed = [3.5, 1.2, 7.8, -2.4, 5.5, 0.0, 9.1, 4.2]
    x_skew = [1.0, 2.0, 3.0, 4.0, 10.0]
    x_norm = rng.standard_normal(500).tolist()

    cor_x = [1.0, 2.0, 3.0, 4.0, 5.0]
    cor_y_pos = [2.0, 4.0, 6.0, 8.0, 10.0]
    cor_y_neg = [10.0, 8.0, 6.0, 4.0, 2.0]
    cor_y_nonlin = [1.0, 4.0, 9.0, 16.0, 25.0]
    x_noise = rng.standard_normal(100).tolist()
    y_noise = (np.array(x_noise) * 0.7 + rng.standard_normal(100) * 0.3).tolist()

    X_mat = rng.standard_normal((50, 3)).tolist()

    data = {}

    data["mean"] = [
        {"x": x_small, "expected": float(np.mean(x_small))},
        {"x": x_mixed, "expected": float(np.mean(x_mixed))},
        {"x": x_norm, "expected": float(np.mean(x_norm))},
    ]

    data["median"] = [
        {"x": x_small, "expected": float(np.median(x_small))},
        {"x": x_med, "expected": float(np.median(x_med))},
        {"x": x_mixed, "expected": float(np.median(x_mixed))},
    ]

    data["mode"] = [
        {"x": [1.0, 2.0, 2.0, 3.0, 3.0, 3.0], "expected": [3.0]},
        {"x": [1.0, 1.0, 2.0, 2.0, 3.0], "expected": [1.0, 2.0]},
        {"x": [5.0, 5.0, 5.0, 1.0, 2.0], "expected": [5.0]},
    ]

    data["trimmed_mean"] = [
        {"x": x_med, "trim": 0.0, "expected": float(st.trim_mean(x_med, 0.0))},
        {"x": x_med, "trim": 0.1, "expected": float(st.trim_mean(x_med, 0.1))},
        {"x": x_med, "trim": 0.2, "expected": float(st.trim_mean(x_med, 0.2))},
        {"x": x_norm, "trim": 0.1, "expected": float(st.trim_mean(x_norm, 0.1))},
    ]

    data["quantile"] = [
        {"x": x_med, "p": 0.0, "expected": float(np.quantile(x_med, 0.0))},
        {"x": x_med, "p": 0.25, "expected": float(np.quantile(x_med, 0.25))},
        {"x": x_med, "p": 0.5, "expected": float(np.quantile(x_med, 0.5))},
        {"x": x_med, "p": 0.75, "expected": float(np.quantile(x_med, 0.75))},
        {"x": x_med, "p": 1.0, "expected": float(np.quantile(x_med, 1.0))},
        {"x": x_norm, "p": 0.1, "expected": float(np.quantile(x_norm, 0.1))},
        {"x": x_norm, "p": 0.9, "expected": float(np.quantile(x_norm, 0.9))},
    ]

    data["percentile"] = [
        {"x": x_med, "p": 25.0, "expected": float(np.percentile(x_med, 25))},
        {"x": x_med, "p": 50.0, "expected": float(np.percentile(x_med, 50))},
        {"x": x_med, "p": 75.0, "expected": float(np.percentile(x_med, 75))},
    ]

    data["var"] = [
        {"x": x_small, "expected": float(np.var(x_small, ddof=1))},
        {"x": x_mixed, "expected": float(np.var(x_mixed, ddof=1))},
        {"x": x_norm, "expected": float(np.var(x_norm, ddof=1))},
    ]

    data["sd"] = [
        {"x": x_small, "expected": float(np.std(x_small, ddof=1))},
        {"x": x_mixed, "expected": float(np.std(x_mixed, ddof=1))},
        {"x": x_norm, "expected": float(np.std(x_norm, ddof=1))},
    ]

    data["iqr"] = [
        {"x": x_med, "expected": float(st.iqr(x_med))},
        {"x": x_norm, "expected": float(st.iqr(x_norm))},
    ]

    data["mad"] = [
        {"x": x_med, "expected": float(st.median_abs_deviation(x_med, scale="normal"))},
        {"x": x_norm, "expected": float(st.median_abs_deviation(x_norm, scale="normal"))},
    ]

    data["cv"] = [
        {"x": x_small, "expected": float(np.std(x_small, ddof=1) / abs(np.mean(x_small)))},
        {"x": x_med, "expected": float(np.std(x_med, ddof=1) / abs(np.mean(x_med)))},
    ]

    data["skewness"] = [
        {"x": x_med, "expected": float(st.skew(x_med, bias=False))},
        {"x": x_skew, "expected": float(st.skew(x_skew, bias=False))},
        {"x": x_norm, "expected": float(st.skew(x_norm, bias=False))},
    ]

    data["kurtosis"] = [
        {"x": x_med + x_med, "expected": float(st.kurtosis(x_med + x_med, bias=False, fisher=True))},
        {"x": x_norm, "expected": float(st.kurtosis(x_norm, bias=False, fisher=True))},
    ]

    data["cor"] = [
        {"x": cor_x, "y": cor_y_pos, "expected": 1.0},
        {"x": cor_x, "y": cor_y_neg, "expected": -1.0},
        {"x": x_noise, "y": y_noise, "expected": float(st.pearsonr(x_noise, y_noise).statistic)},
    ]

    data["spearman"] = [
        {"x": cor_x, "y": cor_y_nonlin, "expected": 1.0},
        {"x": x_noise, "y": y_noise, "expected": float(st.spearmanr(x_noise, y_noise).statistic)},
    ]

    data["kendall"] = [
        {"x": cor_x, "y": cor_y_pos, "expected": 1.0},
        {"x": x_noise[:30], "y": y_noise[:30],
         "expected": float(st.kendalltau(x_noise[:30], y_noise[:30]).statistic)},
    ]

    data["cov"] = [
        {"x": cor_x, "y": cor_y_pos, "expected": float(np.cov(cor_x, cor_y_pos, ddof=1)[0, 1])},
        {"x": x_noise, "y": y_noise, "expected": float(np.cov(x_noise, y_noise, ddof=1)[0, 1])},
    ]

    Xm = np.array(X_mat)
    data["cov_matrix"] = {
        "X": X_mat,
        "expected": np.cov(Xm, rowvar=False, ddof=1).tolist(),
    }

    data["cor_matrix"] = {
        "X": X_mat,
        "expected": np.corrcoef(Xm, rowvar=False).tolist(),
    }

    with open("descriptive.json", "w") as f:
        json.dump(data, f, indent=2)
    n_cases = sum(len(v) if isinstance(v, list) else 1 for v in data.values())
    print(f"  descriptive.json: {n_cases} test groups across {len(data)} functions")


def gen_htest():
    """Generate reference values for Wave 4 parametric hypothesis tests."""
    rng = np.random.default_rng(42)
    data = {}

    # ----- t-tests -----
    x1 = [1.0, 2.0, 3.0, 4.0, 5.0]
    x2 = rng.standard_normal(30).tolist()
    x_norm0 = rng.standard_normal(50).tolist()
    x_norm2 = (np.array(rng.standard_normal(50)) + 2.0).tolist()
    paired_x = [10.0, 12.0, 9.0, 11.0, 13.0, 14.0, 12.0, 11.0]
    paired_y = [9.5, 11.5, 8.5, 10.5, 12.5, 13.0, 11.5, 10.5]

    def ttest1_ref(x, mu0):
        res = st.ttest_1samp(x, mu0)
        return {"statistic": float(res.statistic), "p_value": float(res.pvalue),
                "df": float(len(x) - 1)}

    def ttest2_ref(x, y):
        res = st.ttest_ind(x, y, equal_var=True)
        return {"statistic": float(res.statistic), "p_value": float(res.pvalue),
                "df": float(len(x) + len(y) - 2)}

    def welch_ref(x, y):
        res = st.ttest_ind(x, y, equal_var=False)
        return {"statistic": float(res.statistic), "p_value": float(res.pvalue),
                "df": float(res.df)}

    def ttest_paired_ref(x, y):
        res = st.ttest_rel(x, y)
        return {"statistic": float(res.statistic), "p_value": float(res.pvalue),
                "df": float(len(x) - 1)}

    data["ttest1"] = [
        {"x": x1, "mu0": 3.0, **ttest1_ref(x1, 3.0)},
        {"x": x1, "mu0": 5.0, **ttest1_ref(x1, 5.0)},
        {"x": x2, "mu0": 0.0, **ttest1_ref(x2, 0.0)},
        {"x": x2, "mu0": 0.5, **ttest1_ref(x2, 0.5)},
    ]
    data["ttest2"] = [
        {"x": x_norm0, "y": x_norm2, **ttest2_ref(x_norm0, x_norm2)},
        {"x": x1, "y": [2.0, 3.0, 4.0, 5.0, 6.0], **ttest2_ref(x1, [2.0, 3.0, 4.0, 5.0, 6.0])},
    ]
    data["welch"] = [
        {"x": x_norm0, "y": x_norm2, **welch_ref(x_norm0, x_norm2)},
        {"x": x1, "y": [2.0, 3.0, 4.0, 5.0, 6.0], **welch_ref(x1, [2.0, 3.0, 4.0, 5.0, 6.0])},
        {"x": x_norm0, "y": rng.standard_normal(80).tolist(),
         "_meta": "different n"},
    ]
    # populate the placeholder
    data["welch"][2].update(welch_ref(data["welch"][2]["x"], data["welch"][2]["y"]))
    del data["welch"][2]["_meta"]

    data["ttest_paired"] = [
        {"x": paired_x, "y": paired_y, **ttest_paired_ref(paired_x, paired_y)},
    ]

    # ----- F-test for variance ratio -----
    def ftest_ref(x, y):
        vx = np.var(x, ddof=1)
        vy = np.var(y, ddof=1)
        F = vx / vy
        df1, df2 = len(x) - 1, len(y) - 1
        cdf = st.f.cdf(F, df1, df2)
        pval = 2.0 * min(cdf, 1.0 - cdf)
        return {"statistic": float(F), "df1": df1, "df2": df2, "p_value": float(pval)}

    data["ftest"] = [
        {"x": x_norm0, "y": x_norm2, **ftest_ref(x_norm0, x_norm2)},
        {"x": x1, "y": [10.0, 11.0, 9.0, 12.0, 8.0], **ftest_ref(x1, [10.0, 11.0, 9.0, 12.0, 8.0])},
    ]

    # ----- Chi-squared -----
    def gof_ref(obs, exp):
        # scipy.stats.chisquare needs sum(obs) == sum(exp)
        res = st.chisquare(obs, f_exp=exp)
        return {"statistic": float(res.statistic), "df": len(obs) - 1,
                "p_value": float(res.pvalue)}

    def ind_ref(table):
        chi2, p, dof, _ = st.chi2_contingency(table)
        return {"statistic": float(chi2), "df": int(dof), "p_value": float(p)}

    # scipy.chisquare insists sum(obs) == sum(exp) exactly, so choose totals
    # that divide evenly into k categories.
    obs1 = [10.0, 12.0, 11.0, 9.0, 13.0, 17.0]   # sum 72, k=6, expected=12 each
    exp1 = [12.0] * 6
    obs2 = [27.0, 30.0, 45.0]                    # sum 102, k=3, expected=34 each
    exp2 = [34.0] * 3
    obs3 = [85.0, 95.0, 110.0, 90.0]             # sum 380, k=4, expected=95 each
    exp3 = [95.0] * 4
    data["chisq_gof"] = [
        {"observed": obs1, "expected": exp1, **gof_ref(obs1, exp1)},
        {"observed": obs2, "expected": exp2, **gof_ref(obs2, exp2)},
        {"observed": obs3, "expected": exp3, **gof_ref(obs3, exp3)},
    ]

    table_2x2 = [[10, 15], [20, 25]]
    table_3x3 = [[10, 20, 15], [12, 18, 20], [15, 22, 18]]
    data["chisq_ind"] = [
        {"table": table_2x2, **ind_ref(table_2x2)},
        {"table": table_3x3, **ind_ref(table_3x3)},
    ]

    # ----- ANOVA -----
    def anova1_ref(y, grp):
        groups_split = {}
        for yi, gi in zip(y, grp):
            groups_split.setdefault(gi, []).append(yi)
        res = st.f_oneway(*groups_split.values())
        k = len(groups_split)
        n = len(y)
        return {"statistic": float(res.statistic), "df1": k - 1, "df2": n - k,
                "p_value": float(res.pvalue)}

    anova_y = (10 * [5.0]) + (10 * [10.0]) + (10 * [15.0])
    # add noise so the F-stat is finite
    anova_y = [v + n for v, n in zip(anova_y, rng.standard_normal(30).tolist())]
    anova_grp = (10 * [1]) + (10 * [2]) + (10 * [3])

    anova2_y = rng.standard_normal(60).tolist()
    anova2_grp = (15 * [1]) + (15 * [2]) + (15 * [3]) + (15 * [4])

    data["anova1"] = [
        {"y": anova_y, "grp": anova_grp, **anova1_ref(anova_y, anova_grp)},
        {"y": anova2_y, "grp": anova2_grp, **anova1_ref(anova2_y, anova2_grp)},
    ]

    # ----- Correlation test -----
    def cortest_ref(x, y):
        res = st.pearsonr(x, y)
        n = len(x)
        r = res.statistic
        # t = r * sqrt(n-2) / sqrt(1 - r^2)
        tstat = r * np.sqrt(n - 2) / np.sqrt(1.0 - r * r)
        return {"r": float(r), "statistic": float(tstat),
                "df": n - 2, "p_value": float(res.pvalue)}

    cor_x = list(range(100))
    cor_y = (np.array(cor_x) + rng.standard_normal(100) * 5.0).tolist()
    cor_x2 = rng.standard_normal(50).tolist()
    cor_y2 = rng.standard_normal(50).tolist()

    data["cortest"] = [
        {"x": cor_x, "y": cor_y, **cortest_ref(cor_x, cor_y)},
        {"x": cor_x2, "y": cor_y2, **cortest_ref(cor_x2, cor_y2)},
    ]

    # ----- Proportion test -----
    def proptest_ref(x, n, p0):
        phat = x / n
        se = np.sqrt(p0 * (1 - p0) / n)
        z = (phat - p0) / se
        p = 2.0 * (1.0 - st.norm.cdf(abs(z)))
        return {"statistic": float(z), "p_value": float(p)}

    data["proptest"] = [
        {"x": 60, "n": 100, "p0": 0.5, **proptest_ref(60, 100, 0.5)},
        {"x": 45, "n": 100, "p0": 0.5, **proptest_ref(45, 100, 0.5)},
        {"x": 530, "n": 1000, "p0": 0.5, **proptest_ref(530, 1000, 0.5)},
    ]

    with open("htest.json", "w") as f:
        json.dump(data, f, indent=2)
    n_cases = sum(len(v) for v in data.values())
    print(f"  htest.json: {n_cases} test groups across {len(data)} functions")


def gen_nonparam():
    """Generate reference values for Wave 5 nonparametric tests."""
    rng = np.random.default_rng(42)
    data = {}

    # Test samples
    x_norm = rng.standard_normal(50).tolist()
    y_norm = rng.standard_normal(50).tolist()
    y_shift = (np.array(rng.standard_normal(50)) + 1.5).tolist()
    x_small = [3.0, 5.0, 7.0, 9.0, 11.0]
    y_small = [4.0, 6.0, 8.0, 10.0, 12.0]
    x_ties = [1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 5.0]
    y_ties = [2.0, 2.0, 3.0, 4.0, 4.0, 5.0, 6.0, 7.0]

    # Paired samples for Wilcoxon
    paired_x = rng.standard_normal(30).tolist()
    paired_y = (np.array(paired_x) + 0.5 + rng.standard_normal(30) * 0.2).tolist()
    paired_x_ties = [10.0, 12.0, 9.0, 11.0, 13.0, 14.0, 12.0, 11.0, 10.0, 12.0]
    paired_y_ties = [9.5, 11.5, 9.0, 10.5, 12.5, 13.0, 11.5, 10.5, 9.5, 11.5]

    # ----- KS test -----
    # NOTE: scipy.ks_2samp(method='asymp') actually uses kstwo.sf (the FINITE-n
    # Kolmogorov-Smirnov distribution), not the limiting Kolmogorov SF that the
    # plan calls for. We reference scipy.special.kolmogorov (the limiting form,
    # SF(y) = 2 sum (-1)^(k-1) exp(-2 k^2 y^2)) directly so the q impl can match.
    def ks_ref(x, y):
        nx, ny = len(x), len(y)
        x_arr, y_arr = np.array(x), np.array(y)
        # Match the q impl's D computation on the pooled support
        pooled = np.unique(np.concatenate([x_arr, y_arr]))
        fx = np.array([(x_arr <= t).mean() for t in pooled])
        fy = np.array([(y_arr <= t).mean() for t in pooled])
        d = float(np.max(np.abs(fx - fy)))
        en = np.sqrt(nx * ny / (nx + ny))
        pval = float(sp.kolmogorov(en * d))
        return {"statistic": d, "p_value": max(0.0, min(1.0, pval))}

    data["ks"] = [
        {"x": x_norm, "y": y_norm, **ks_ref(x_norm, y_norm)},
        {"x": x_norm, "y": y_shift, **ks_ref(x_norm, y_shift)},
        {"x": x_small, "y": y_small, **ks_ref(x_small, y_small)},
    ]

    # ----- Mann-Whitney U -----
    def mw_ref(x, y):
        res = st.mannwhitneyu(x, y, alternative="two-sided",
                              method="asymptotic", use_continuity=True)
        return {"statistic": float(res.statistic), "p_value": float(res.pvalue)}

    data["mannwhitney"] = [
        {"x": x_norm, "y": y_norm, **mw_ref(x_norm, y_norm)},
        {"x": x_norm, "y": y_shift, **mw_ref(x_norm, y_shift)},
        {"x": x_small, "y": y_small, **mw_ref(x_small, y_small)},
        {"x": x_ties, "y": y_ties, **mw_ref(x_ties, y_ties)},
    ]

    # ----- Wilcoxon signed-rank -----
    def wcx_ref(x, y):
        # scipy returns min(W+, W-) as statistic in two-sided mode; we will
        # match by storing both W+ and the scipy statistic.
        res = st.wilcoxon(x, y, alternative="two-sided", method="approx",
                          correction=True, zero_method="wilcox")
        d = np.array(x) - np.array(y)
        d = d[d != 0]
        ranks_abs = st.rankdata(np.abs(d))
        w_plus = float(np.sum(ranks_abs[d > 0]))
        return {"statistic": float(res.statistic), "w_plus": w_plus,
                "p_value": float(res.pvalue)}

    data["wilcoxon"] = [
        {"x": paired_x, "y": paired_y, **wcx_ref(paired_x, paired_y)},
        {"x": paired_x_ties, "y": paired_y_ties, **wcx_ref(paired_x_ties, paired_y_ties)},
    ]

    # ----- Shapiro-Wilk -----
    def shapiro_ref(x):
        res = st.shapiro(x)
        return {"statistic": float(res.statistic), "p_value": float(res.pvalue)}

    sw_norm15 = rng.standard_normal(15).tolist()
    sw_norm30 = rng.standard_normal(30).tolist()
    sw_norm100 = rng.standard_normal(100).tolist()
    sw_unif50 = rng.uniform(-1, 1, 50).tolist()
    sw_exp50 = rng.exponential(1.0, 50).tolist()

    data["shapiro"] = [
        {"x": sw_norm15, **shapiro_ref(sw_norm15)},
        {"x": sw_norm30, **shapiro_ref(sw_norm30)},
        {"x": sw_norm100, **shapiro_ref(sw_norm100)},
        {"x": sw_unif50, **shapiro_ref(sw_unif50)},
        {"x": sw_exp50, **shapiro_ref(sw_exp50)},
    ]

    # Add small-sample Shapiro-Wilk cases (n in [4, 11]) for Wave 8 Gap B2
    sw_small_samples = {
        4: rng.standard_normal(4).tolist(),
        5: rng.standard_normal(5).tolist(),
        6: rng.standard_normal(6).tolist(),
        7: rng.standard_normal(7).tolist(),
        8: rng.standard_normal(8).tolist(),
        9: rng.standard_normal(9).tolist(),
        10: rng.standard_normal(10).tolist(),
        11: rng.standard_normal(11).tolist(),
    }
    data["shapiro_small"] = [
        {"n": n, "x": x, **shapiro_ref(x)}
        for n, x in sw_small_samples.items()
    ]
    # Also add edge cases: uniform and exponential for small n
    sw_unif8 = rng.uniform(-1, 1, 8).tolist()
    sw_exp10 = rng.exponential(1.0, 10).tolist()
    data["shapiro_small"].extend([
        {"n": 8, "x": sw_unif8, **shapiro_ref(sw_unif8)},
        {"n": 10, "x": sw_exp10, **shapiro_ref(sw_exp10)},
    ])

    # ----- Jarque-Bera -----
    def jb_ref(x):
        res = st.jarque_bera(x)
        return {"statistic": float(res.statistic), "df": 2,
                "p_value": float(res.pvalue)}

    data["jarque_bera"] = [
        {"x": rng.standard_normal(100).tolist(),
         **jb_ref(rng.standard_normal(100).tolist())},
        # NOTE: above generates fresh samples for each call; instead use the
        # same sample for x and ref:
    ]
    jb_samples = [
        rng.standard_normal(100).tolist(),
        rng.exponential(1.0, 100).tolist(),
        rng.standard_normal(500).tolist(),
        list(range(20)),
    ]
    data["jarque_bera"] = [{"x": s, **jb_ref(s)} for s in jb_samples]

    with open("nonparam.json", "w") as f:
        json.dump(data, f, indent=2)
    n_cases = sum(len(v) for v in data.values())
    print(f"  nonparam.json: {n_cases} test groups across {len(data)} functions")


def gen_diagnostics():
    """Generate reference values for Wave 6 regression diagnostics."""
    import statsmodels.api as sm
    from statsmodels.stats.diagnostic import het_breuschpagan, het_white
    from statsmodels.stats.stattools import durbin_watson
    from statsmodels.stats.outliers_influence import variance_inflation_factor

    rng = np.random.default_rng(42)
    n = 30
    p_nonconst = 3
    # Build design matrix with intercept + 3 predictors, then deliberately
    # introduce one high-leverage point and one outlier to exercise the
    # influence diagnostics.
    X_raw = rng.standard_normal((n, p_nonconst))
    X_raw[5, :] *= 4.0  # high leverage row
    X = np.hstack([np.ones((n, 1)), X_raw])
    true_beta = np.array([1.5, 0.7, -0.4, 1.1])
    y = X @ true_beta + rng.standard_normal(n) * 0.5
    y[10] += 4.0  # outlier on the response

    model = sm.OLS(y, X).fit()
    infl = model.get_influence()

    data = {}
    data["dataset"] = {
        "X": X.tolist(),
        "y": y.tolist(),
        "beta": model.params.tolist(),
        "residuals": model.resid.tolist(),
        "fitted": model.fittedvalues.tolist(),
        "leverage": infl.hat_matrix_diag.tolist(),
        "n": int(n),
        "p": int(X.shape[1]),
        "rss": float(model.ssr),
        "tss": float(model.centered_tss),
    }

    # VIF: statsmodels returns VIF per column including intercept; we skip
    # the intercept column (R `car::vif` convention, matches the plan).
    data["vif"] = [float(variance_inflation_factor(X, i))
                   for i in range(1, X.shape[1])]

    data["cooks_distance"] = infl.cooks_distance[0].tolist()
    data["dffits"] = infl.dffits[0].tolist()
    data["dfbetas"] = infl.dfbetas.tolist()

    data["durbin_watson"] = float(durbin_watson(model.resid))

    # Durbin-Watson p-value via Pan's beta approximation (eigenvalue-based)
    # statsmodels doesn't compute this, so implement it here
    def durbin_watson_pvalue(dw_stat, X, n, p):
        # Construct A (differencing matrix)
        A = np.zeros((n, n))
        for i in range(1, n):
            A[i, i] = 1
            A[i, i-1] = -1
        # Compute M = I - H
        H = X @ np.linalg.inv(X.T @ X) @ X.T
        M = np.eye(n) - H
        # Eigenvalues of M @ A @ M
        MAM = M @ A @ M
        eigvals = np.linalg.eigvalsh(MAM)  # symmetric, so eigvalsh
        # Mean and variance
        mu = eigvals.sum() / (n - p)
        var = 2 * (np.sum(eigvals**2) - eigvals.sum()**2 / (n - p)) / (n - p)**2
        # Guard against zero variance
        if var < 1e-12:
            var = 1e-12
        # Beta parameters
        a = mu * (mu * (4 - mu) / var - 1) / 2
        b = (4 - mu) * a / mu if mu > 1e-12 else a
        # p-value
        x = dw_stat / 4.0
        from scipy.stats import beta
        pval = beta.cdf(x, a, b)
        return float(pval)

    data["durbin_watson_pvalue"] = durbin_watson_pvalue(
        data["durbin_watson"], X, n, X.shape[1]
    )

    bp = het_breuschpagan(model.resid, X)
    # statsmodels het_breuschpagan returns (lm_stat, lm_pvalue, fvalue, f_pvalue)
    data["breusch_pagan"] = {"statistic": float(bp[0]), "df": X.shape[1] - 1,
                             "p_value": float(bp[1])}

    wh = het_white(model.resid, X)
    data["white_test"] = {"statistic": float(wh[0]),
                          "p_value": float(wh[1])}

    data["rsquared"] = {"rsq": float(model.rsquared),
                        "adj_rsq": float(model.rsquared_adj)}
    data["aic"] = float(model.aic)
    data["bic"] = float(model.bic)

    # Jarque-Bera on residuals via the same Wave 5 scipy.stats reference,
    # except we use biased moments (matches scipy.stats.jarque_bera).
    jb_res = st.jarque_bera(model.resid)
    data["jarque_bera"] = {"statistic": float(jb_res.statistic),
                           "df": 2, "p_value": float(jb_res.pvalue)}

    with open("diagnostics.json", "w") as f:
        json.dump(data, f, indent=2)
    print(f"  diagnostics.json: 1 dataset + {len(data) - 1} diagnostic groups")


def gen_dist_extra():
    """Generate reference values for Wave 7 additional distributions."""
    data = {}

    # ---------- Beta ----------
    beta_params = [(0.5, 0.5), (2.0, 5.0), (5.0, 2.0), (2.0, 2.0)]
    x_beta = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95]
    p_beta = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
    beta_groups = {}
    for a, b in beta_params:
        key = f"a{a}_b{b}"
        beta_groups[key] = {
            "alpha": a, "beta": b,
            "dbeta": {"x": x_beta, "y": [float(st.beta.pdf(x, a, b)) for x in x_beta]},
            "pbeta": {"x": x_beta, "y": [float(st.beta.cdf(x, a, b)) for x in x_beta]},
            "qbeta": {"p": p_beta, "y": [float(st.beta.ppf(p, a, b)) for p in p_beta]},
        }
    data["beta"] = beta_groups

    # ---------- Gamma (shape-rate; scipy uses scale = 1/rate) ----------
    gamma_params = [(1.0, 1.0), (2.0, 0.5), (5.0, 1.0), (0.5, 2.0)]
    x_gamma = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    p_gamma = [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
    gamma_groups = {}
    for a, b in gamma_params:
        key = f"a{a}_b{b}"
        scale = 1.0 / b
        gamma_groups[key] = {
            "alpha": a, "beta": b,
            "dgamma": {"x": x_gamma, "y": [float(st.gamma.pdf(x, a, scale=scale)) for x in x_gamma]},
            "pgamma": {"x": x_gamma, "y": [float(st.gamma.cdf(x, a, scale=scale)) for x in x_gamma]},
            "qgamma": {"p": p_gamma, "y": [float(st.gamma.ppf(p, a, scale=scale)) for p in p_gamma]},
        }
    data["gamma"] = gamma_groups

    # ---------- Binomial ----------
    binom_params = [(10, 0.5), (20, 0.3), (50, 0.1), (100, 0.5)]
    binom_groups = {}
    for n, p in binom_params:
        ks = list(range(0, n + 1, max(1, n // 10)))
        ps = [0.01, 0.1, 0.5, 0.9, 0.99]
        binom_groups[f"n{n}_p{p}"] = {
            "n": n, "p": p,
            "dbinom": {"k": ks, "y": [float(st.binom.pmf(k, n, p)) for k in ks]},
            "pbinom": {"k": ks, "y": [float(st.binom.cdf(k, n, p)) for k in ks]},
            "qbinom": {"prob": ps, "y": [int(st.binom.ppf(p2, n, p)) for p2 in ps]},
        }
    data["binom"] = binom_groups

    # ---------- Poisson ----------
    pois_lams = [0.5, 1.0, 3.0, 10.0, 50.0]
    pois_groups = {}
    for lam in pois_lams:
        ks = list(range(0, max(10, int(3 * lam) + 1)))
        ps = [0.01, 0.1, 0.5, 0.9, 0.99]
        pois_groups[f"lam{lam}"] = {
            "lambda": lam,
            "dpois": {"k": ks, "y": [float(st.poisson.pmf(k, lam)) for k in ks]},
            "ppois": {"k": ks, "y": [float(st.poisson.cdf(k, lam)) for k in ks]},
            "qpois": {"p": ps, "y": [int(st.poisson.ppf(p, lam)) for p in ps]},
        }
    data["poisson"] = pois_groups

    # ---------- Exponential ----------
    exp_rates = [0.5, 1.0, 2.0, 5.0]
    x_exp = [0.1, 0.5, 1.0, 2.0, 5.0]
    p_exp = [0.01, 0.1, 0.5, 0.9, 0.99]
    exp_groups = {}
    for rate in exp_rates:
        scale = 1.0 / rate
        exp_groups[f"rate{rate}"] = {
            "rate": rate,
            "dexp": {"x": x_exp, "y": [float(st.expon.pdf(x, scale=scale)) for x in x_exp]},
            "pexp": {"x": x_exp, "y": [float(st.expon.cdf(x, scale=scale)) for x in x_exp]},
            "qexp": {"p": p_exp, "y": [float(st.expon.ppf(p, scale=scale)) for p in p_exp]},
        }
    data["expon"] = exp_groups

    with open("dist_extra.json", "w") as f:
        json.dump(data, f, indent=2)
    n_dists = len(data)
    n_groups = sum(len(g) for g in data.values())
    print(f"  dist_extra.json: {n_groups} groups across {n_dists} distributions")


if __name__ == "__main__":
    print("Generating qstats reference values...")
    gen_special()
    gen_distributions()
    gen_linalg()
    gen_descriptive()
    gen_htest()
    gen_nonparam()
    gen_diagnostics()
    gen_dist_extra()
    print("Done. All reference JSON files generated.")
