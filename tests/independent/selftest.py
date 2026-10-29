"""Synthetic self-checks of my own CR2 code (no project data): vectorised cr2_test vs a brute-force full-matrix CR2
(Bell-McCaffrey, explicit (I - H_cc)^-1/2), and vs statsmodels HC2 when every cluster is a singleton."""
from __future__ import annotations

import numpy as np

from .confounder import cr2_test, orth, resid


def brute(Z, d, y, cl):
    X = np.column_stack([Z, d])
    XtXi = np.linalg.pinv(X.T @ X)
    b = XtXi @ X.T @ y
    e = y - X @ b
    H = X @ XtXi @ X.T
    meat = np.zeros((X.shape[1], X.shape[1]))
    for c in np.unique(cl):
        i = np.where(cl == c)[0]
        w, V = np.linalg.eigh(np.eye(len(i)) - H[np.ix_(i, i)])
        A = (V / np.sqrt(w)) @ V.T
        u = X[i].T @ A @ e[i]
        meat += np.outer(u, u)
    Vb = XtXi @ meat @ XtXi
    return b[-1], np.sqrt(Vb[-1, -1])


def run(seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for rep in range(20):
        n = 200
        Z = np.column_stack([np.ones(n), rng.normal(size=(n, 4)), (rng.random(n) < 0.3).astype(float)])
        D = np.column_stack([(rng.random(n) < 0.05 + 0.05 * k).astype(float) for k in range(3)] + [rng.normal(size=n)])
        y = rng.normal(size=n) * (1 + D[:, 0])
        cl = np.arange(n).astype(str)
        if rep % 2:
            cl[:30] = np.repeat(np.arange(10), 3).astype(str)  # some 3-line clusters
            cl[30:40] = np.repeat(np.arange(100, 105), 2).astype(str)
        Q = orth(Z)
        b, s, _ = cr2_test(Q, resid(Q, y[:, None])[:, 0], resid(Q, D), cl)
        for j in range(D.shape[1]):
            bb, ss = brute(Z, D[:, j], y, cl)
            out.append(max(abs(b[j] - bb) / abs(bb), abs(s[j] - ss) / ss))
        if rep % 2 == 0:
            import statsmodels.api as sm
            for j in range(D.shape[1]):
                r = sm.OLS(y, np.column_stack([Z, D[:, j]])).fit(cov_type="HC2")
                out.append(abs(np.sqrt(r.cov_params()[-1, -1]) - s[j]) / s[j])
    return float(max(out))


if __name__ == "__main__":
    print("max rel diff", run())
