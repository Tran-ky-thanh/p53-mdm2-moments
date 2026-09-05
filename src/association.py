# -*- coding: utf-8 -*-
"""
Utilities: (1) add scRNA-seq-style NOISE, (2) DEPENDENCE measures between two variables.

Noise:
  - dropout/capture: each mRNA molecule is "captured" with probability beta (binomial thinning)
    -> models the low capture efficiency & dropout of scRNA-seq (a caveat noted in the paper).
  - (optional) additive Gaussian noise on log1p.

Dependence measures between X (p53 mRNA) and Y (MDM2 mRNA) across the cell population:
  - Pearson r (linear), Spearman rho (monotone)
  - Mutual Information (kNN/Kraskov, via sklearn) -- captures nonlinear dependence
  - normalized HSIC (RBF kernel, median heuristic) -- captures general nonlinear dependence
"""
import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.feature_selection import mutual_info_regression
try:
    import dcor as _dcor
except ImportError:
    _dcor = None


def add_dropout(counts, beta, rng):
    """Binomial thinning: keep each molecule with probability beta (0<beta<=1)."""
    counts = np.asarray(counts)
    if beta >= 1.0:
        return counts.astype(float)
    return rng.binomial(counts.astype(int), beta).astype(float)


def splatter_noise(C, rng, lib_scale=0.2, bcv=0.5, dropout_mid=1.0, dropout_shape=-1.0):
    """
    SPLATTER-style scRNA-seq technical-noise model (Zappia et al. 2017, Genome Biol 18:174) --
    the most widely used scRNA-seq simulator. Applied to "true" counts C (N_cells, N_genes) from SSA.

    Components:
      1) per-cell library size ~ log-normal (variable sequencing depth).
      2) Negative-Binomial overdispersion via BCV (gamma-Poisson) -> mean-variance trend.
      3) EXPRESSION-DEPENDENT DROPOUT: zero probability is a logistic function of log(mean);
         lowly-expressed genes drop out more (a hallmark of scRNA-seq).
    """
    C = np.asarray(C, float)
    n_cells = C.shape[0]
    s = rng.lognormal(mean=-0.5 * lib_scale**2, sigma=lib_scale, size=(n_cells, 1))  # mean~1
    base = C * s
    if bcv and bcv > 0:
        shape = 1.0 / (bcv**2)
        rate = np.where(base > 0, rng.gamma(shape=shape, scale=np.maximum(base, 1e-12) / shape), 0.0)
    else:
        rate = base
    counts = rng.poisson(rate).astype(float)
    logmean = np.log(base + 1e-8)
    pi = 1.0 / (1.0 + np.exp(-dropout_shape * (logmean - dropout_mid)))  # dropout probability
    keep = rng.random(counts.shape) > pi
    return counts * keep


def add_gaussian_log(counts, sigma, rng):
    """Gaussian noise on log1p, returned on the count scale (>=0)."""
    z = np.log1p(counts) + rng.normal(0, sigma, size=counts.shape)
    return np.clip(np.expm1(z), 0, None)


# ------------------------- dependence measures -------------------------
def pearson(x, y):
    if np.std(x) == 0 or np.std(y) == 0:
        return 0.0
    return float(pearsonr(x, y)[0])


def spearman(x, y):
    if np.std(x) == 0 or np.std(y) == 0:
        return 0.0
    return float(spearmanr(x, y)[0])


def mutual_info(x, y, seed=0):
    """MI(X;Y) [nats] via sklearn kNN (Kraskov) estimator."""
    return float(mutual_info_regression(x.reshape(-1, 1), y, random_state=seed)[0])


def _rbf_gram(x):
    x = np.asarray(x, float).reshape(-1, 1)
    d2 = (x - x.T) ** 2
    med = np.median(d2[d2 > 0]) if np.any(d2 > 0) else 1.0
    sigma2 = 0.5 * med                      # median heuristic
    return np.exp(-d2 / (2.0 * sigma2 + 1e-12))


def hsic_normalized(x, y):
    """
    Normalized HSIC (nHSIC) in [0,1]: 0 = independent, larger = more dependent.
    nHSIC = HSIC(K,L) / sqrt(HSIC(K,K)*HSIC(L,L)) with an RBF kernel (median heuristic).
    """
    n = len(x)
    K = _rbf_gram(x); L = _rbf_gram(y)
    H = np.eye(n) - 1.0 / n
    Kc = H @ K @ H; Lc = H @ L @ H
    hxy = np.sum(Kc * Lc)
    hxx = np.sum(Kc * Kc); hyy = np.sum(Lc * Lc)
    if hxx <= 0 or hyy <= 0:
        return 0.0
    return float(hxy / np.sqrt(hxx * hyy))


def distance_corr(x, y):
    """Distance correlation (Szekely) via the dcor library: 0 <=> independent; captures nonlinear dependence."""
    if _dcor is None:
        return float('nan')
    return float(_dcor.distance_correlation(np.asarray(x, float), np.asarray(y, float)))


def all_measures(x, y, n_kernel=1200, seed=0):
    """Return dict {pearson, spearman, MI, HSIC, dCor}. HSIC/dCor are subsampled for speed."""
    rng = np.random.default_rng(seed)
    out = dict(pearson=pearson(x, y), spearman=spearman(x, y), MI=mutual_info(x, y, seed))
    if len(x) > n_kernel:
        idx = rng.choice(len(x), n_kernel, replace=False)
        xh, yh = x[idx], y[idx]
    else:
        xh, yh = x, y
    out['HSIC'] = hsic_normalized(xh, yh)
    out['dCor'] = distance_corr(xh, yh)
    return out


# ------------------- permutation test -------------------
def _double_center_dist(x):
    x = np.asarray(x, float).reshape(-1, 1)
    D = np.abs(x - x.T)
    return D - D.mean(0, keepdims=True) - D.mean(1, keepdims=True) + D.mean()


def _center_gram(x):
    K = _rbf_gram(x); n = len(x); H = np.eye(n) - 1.0 / n
    return H @ K @ H


def permutation_z(x, y, method='dCor', B=200, n_sub=800, seed=0):
    """
    Permutation z-score & p-value for the x-y relationship. Uses an index-permutation trick on the
    pre-centered matrices (fast) for HSIC/dCor; direct recomputation for pearson/spearman/MI.
    Returns (observed_stat, z, p_value).
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) > n_sub:
        idx = rng.choice(len(x), n_sub, replace=False); x = x[idx]; y = y[idx]
    n = len(x)
    if method in ('HSIC', 'dCor'):
        if method == 'HSIC':
            Kc = _center_gram(x); Lc = _center_gram(y)
            denom = np.sqrt(np.sum(Kc * Kc) * np.sum(Lc * Lc))
            obs = np.sum(Kc * Lc) / denom
            null = np.empty(B)
            for b in range(B):
                p = rng.permutation(n)
                null[b] = np.sum(Kc * Lc[np.ix_(p, p)]) / denom
        else:
            A = _double_center_dist(x); Bm = _double_center_dist(y)
            dvar = np.sqrt(np.mean(A * A) * np.mean(Bm * Bm))
            obs = np.sqrt(max(np.mean(A * Bm), 0) / dvar) if dvar > 0 else 0.0
            null = np.empty(B)
            for b in range(B):
                p = rng.permutation(n)
                cov = np.mean(A * Bm[np.ix_(p, p)])
                null[b] = np.sqrt(max(cov, 0) / dvar) if dvar > 0 else 0.0
    else:
        func = {'pearson': pearson, 'spearman': spearman,
                'MI': lambda a, b: mutual_info(a, b, seed)}[method]
        obs = func(x, y)
        null = np.array([func(x, rng.permutation(y)) for _ in range(B)])
    mu, sd = null.mean(), null.std() + 1e-12
    z = (obs - mu) / sd
    p = (np.sum(np.abs(null - mu) >= np.abs(obs - mu)) + 1) / (B + 1)
    return float(obs), float(z), float(p)
