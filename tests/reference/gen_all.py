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

    with open("distributions.json", "w") as f:
        json.dump(data, f, indent=2)
    print(f"  distributions.json: normal (dnorm/pnorm/qnorm)")


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


if __name__ == "__main__":
    print("Generating qstats reference values...")
    gen_special()
    gen_distributions()
    gen_linalg()
    print("Done. All reference JSON files generated.")
