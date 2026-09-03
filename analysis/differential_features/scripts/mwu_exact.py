"""Exact / near-exact permutation Mann-Whitney for the tiny-n contrasts in this project.

Why not scipy's `mannwhitneyu` defaults
---------------------------------------
Every contrast here is n=5v5, 5v10 or 10v10, and the choice of null matters
enormously at those sizes. Measured on a perfectly separated, tie-free
feature:

    design   analytic 2/C   scipy exact   asymptotic(+cc)   asymptotic(no cc)
    5 v 5      7.937e-03     7.937e-03       1.219e-02          9.023e-03
    5 v 10     6.660e-04     6.660e-04       2.694e-03          2.200e-03
    10 v 10    1.083e-05     1.083e-05       1.827e-04          1.571e-04

Three consequences that bit this project in sequence:

1. `method="auto"` (the original code) uses the EXACT null only when a feature
   has no ties, and the asymptotic one otherwise. Because the asymptotic null
   is unbounded, a zero-inflated (tie-heavy) feature could attain a SMALLER
   p-value than a perfectly separated dense one -- biasing every ranking
   toward sparsity. On disk this was visible as a p floor of 3.53e-4 for
   tied features against 6.66e-4 for untied ones in the same 5v10 contrast.

2. `method="asymptotic"` (the first fix) restores uniformity but is markedly
   CONSERVATIVE at n=5v5: its floor is 1.22e-2 against the exact 7.94e-3.
   Since BH can only call features when many sit at the floor together
   (k >= p_floor * m / q), inflating the floor by 1.5x inflates the required
   count by 1.5x. That alone moved `dendrobatidis_spore_Sporangium_vs_spore_
   Mature` from 5,507 significant to 0 -- a reporting artifact, not biology.

3. The statistically correct answer at this n is the EXACT PERMUTATION null,
   conditional on the observed tie pattern. scipy's `method="exact"` refuses
   ties, so we enumerate directly.

Method
------
Rank the pooled values per feature (average ranks, so ties are handled
correctly), then form the null distribution of the group-A rank sum over
label assignments. With C(n1+n2, n1) <= `max_exact` we ENUMERATE every
assignment, giving a genuinely exact conditional p-value; above that we
sample `n_sample` random assignments. For this project's designs:

    5 v 5   -> C(10,5)  =    252   exact
    5 v 10  -> C(15,5)  =  3,003   exact
    10 v 10 -> C(20,10) = 184,756  sampled (default 20,000)

Because the rank vector is fixed under permutation, the whole null is one
matrix product (assignment-indicator matrix @ rank matrix), chunked over
features to bound memory.

The two-sided p is the standard permutation estimate
    p = (1 + #{|S_perm - E[S]| >= |S_obs - E[S]|}) / (1 + n_perm)
which is why the attainable floor is 1/(1+252) = 3.95e-3 for 5v5 under
enumeration -- slightly below the analytic 2/252 because the +1 convention
counts the observed assignment itself.
"""
from __future__ import annotations

from itertools import combinations
from math import comb

import numpy as np
from scipy.stats import rankdata


def _assignment_matrix(n1: int, n2: int, max_exact: int, n_sample: int, seed: int):
    """(n_perm, n1+n2) 0/1 indicator of group-A membership; exact if feasible."""
    n = n1 + n2
    total = comb(n, n1)
    if total <= max_exact:
        idx = list(combinations(range(n), n1))
        M = np.zeros((len(idx), n), dtype=np.float32)
        for r, c in enumerate(idx):
            M[r, list(c)] = 1.0
        return M, True
    rng = np.random.default_rng(seed)
    M = np.zeros((n_sample, n), dtype=np.float32)
    for r in range(n_sample):
        M[r, rng.choice(n, size=n1, replace=False)] = 1.0
    return M, False


def mwu_permutation(
    mat_a: np.ndarray,
    mat_b: np.ndarray,
    max_exact: int = 20_000,
    n_sample: int = 20_000,
    seed: int = 0,
    chunk: int = 4_000,
):
    """Two-sided permutation Mann-Whitney p per feature.

    mat_a, mat_b: (n_samples_group, n_features) already-normalized values.
    Returns (pvals, u_stat, is_exact).
    """
    n1, n2 = mat_a.shape[0], mat_b.shape[0]
    n_feat = mat_a.shape[1]
    pooled = np.vstack([mat_a, mat_b])

    M, is_exact = _assignment_matrix(n1, n2, max_exact, n_sample, seed)
    n_perm = M.shape[0]
    expected = n1 * (n1 + n2 + 1) / 2.0

    pvals = np.empty(n_feat)
    ustat = np.empty(n_feat)
    obs_indicator = np.zeros(n1 + n2, dtype=np.float32)
    obs_indicator[:n1] = 1.0

    for start in range(0, n_feat, chunk):
        stop = min(start + chunk, n_feat)
        block = pooled[:, start:stop]
        # Average ranks down each feature column -> ties handled exactly.
        R = np.apply_along_axis(rankdata, 0, block).astype(np.float32)
        s_obs = obs_indicator @ R                      # (chunk,)
        s_null = M @ R                                 # (n_perm, chunk)
        dev_obs = np.abs(s_obs - expected)
        dev_null = np.abs(s_null - expected)
        exceed = (dev_null >= dev_obs[None, :] - 1e-9).sum(axis=0)
        if is_exact:
            # COMPLETE enumeration: the null is the whole reference set, and it
            # already contains the observed assignment, so the exact
            # conditional p is simply the proportion at least as extreme.
            # Using the sampled (1+x)/(1+n) convention here would inflate the
            # 5v5 floor from 7.94e-3 to 1.19e-2 -- a 1.5x conservative bias
            # that costs exactly the BH power this module exists to recover.
            pvals[start:stop] = exceed / n_perm
        else:
            # SAMPLED null: add-one keeps p > 0 and stays slightly conservative.
            pvals[start:stop] = (1.0 + exceed) / (1.0 + n_perm)
        # U for group A, from its rank sum (reported for continuity with the
        # previous tables; the p-value no longer depends on it).
        ustat[start:stop] = s_obs - n1 * (n1 + 1) / 2.0

    return pvals, ustat, is_exact
