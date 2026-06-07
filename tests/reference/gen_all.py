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

    # Log-gamma
    x_lgamma = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 10.0, 50.0, 100.0, 0.01, 0.001]
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

    data["chisq"] = {}
    for df in df_chisq:
        data["chisq"][f"df_{int(df)}"] = {
            "df": df,
            "dchisq": {"x": x_chisq, "y": [float(st.chi2.pdf(x, df)) for x in x_chisq]},
            "pchisq": {"x": x_chisq, "y": [float(st.chi2.cdf(x, df)) for x in x_chisq]},
            "qchisq": {"p": p_chisq, "y": [float(st.chi2.ppf(p, df)) for p in p_chisq]}
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
            "qt": {"p": p_t, "y": [float(st.t.ppf(p, df)) for p in p_t]}
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


if __name__ == "__main__":
    print("Generating qstats reference values...")
    gen_special()
    gen_distributions()
    gen_linalg()
    gen_descriptive()
    print("Done. All reference JSON files generated.")
