# -*- coding: utf-8 -*-
"""
Utilities that turn SYMMETRIC measures (Cor/MI/HSIC/dCor) into DIRECTIONAL ones and separate
REGULATION from CORRELATION. Two principles:

  (1) DIRECTION requires TIME (temporal precedence):
      - lagged cross-correlation / lagged dCor: measured between X(t) and Y(t+tau).
      - Granger causality (1-step): does X(t) improve prediction of Y(t+1) beyond Y(t)?
      - Transfer entropy TE(X->Y) = I(Y_{t+1}; X_t | Y_t) -- the DIRECTIONAL version of MI.
      Asymmetry between directions => who drives whom.

  (2) REGULATION vs CORRELATION requires CONDITIONING:
      - detrending (removing the per-time mean) => removes the shared time trend.
      - partial correlation / conditional MI: condition on a third variable (a common driver)
        or on the target's own past => remove indirect links / common-driver effects.

SSA data provides PER-CELL trajectories, so lagged statistics can be computed (which real
scRNA-seq snapshots do NOT have -- there one must use moment dynamics/MBI or RNA velocity).
"""
import numpy as np
try:
    import dcor as _dcor
except ImportError:
    _dcor = None


def detrend_by_time(X):
    """Remove the population mean at each time -> leaves fluctuations around the shared trend."""
    return X - X.mean(axis=0, keepdims=True)


def lagged_xcorr(X, Y, lag):
    """corr(X(t), Y(t+lag)) pooled over cells & time. lag>0: X leads Y."""
    if lag >= 0:
        x = X[:, :X.shape[1]-lag].ravel(); y = Y[:, lag:].ravel()
    else:
        x = X[:, -lag:].ravel(); y = Y[:, :Y.shape[1]+lag].ravel()
    if x.std() == 0 or y.std() == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def lagged_dcor(X, Y, lag, n_sub=1500, seed=0):
    """Lagged (nonlinear) distance correlation."""
    if _dcor is None:
        return float('nan')
    if lag >= 0:
        x = X[:, :X.shape[1]-lag].ravel(); y = Y[:, lag:].ravel()
    else:
        x = X[:, -lag:].ravel(); y = Y[:, :Y.shape[1]+lag].ravel()
    rng = np.random.default_rng(seed)
    if len(x) > n_sub:
        idx = rng.choice(len(x), n_sub, replace=False); x = x[idx]; y = y[idx]
    return float(_dcor.distance_correlation(x, y))


def granger_1step(src, tgt):
    """1-step linear Granger: log(var_reduced/var_full). >0 => src helps predict tgt."""
    from numpy.linalg import lstsq
    Xt = tgt[:, :-1].ravel(); Xt1 = tgt[:, 1:].ravel(); St = src[:, :-1].ravel()
    base = np.c_[np.ones_like(Xt), Xt]
    full = np.c_[base, St]
    r0 = Xt1 - base.dot(lstsq(base, Xt1, rcond=None)[0])
    r1 = Xt1 - full.dot(lstsq(full, Xt1, rcond=None)[0])
    return float(np.log(max(r0.var(), 1e-12) / max(r1.var(), 1e-12)))


def transfer_entropy(src, tgt, bins=6):
    """
    Transfer entropy TE(src->tgt) = I(tgt_{t+1}; src_t | tgt_t) [nats] (histogram estimator).
    = sum p(y1,y0,x0) log[ p(y1|y0,x0) / p(y1|y0) ].
    """
    y1 = tgt[:, 1:].ravel(); y0 = tgt[:, :-1].ravel(); x0 = src[:, :-1].ravel()
    def binidx(v):
        edges = np.quantile(v, np.linspace(0, 1, bins + 1))
        edges[-1] += 1e-9
        return np.clip(np.digitize(v, edges[1:-1]), 0, bins - 1)
    a = binidx(y1); b = binidx(y0); c = binidx(x0)
    # count joint p(y1,y0,x0)
    P = np.zeros((bins, bins, bins))
    for i in range(len(a)):
        P[a[i], b[i], c[i]] += 1
    P /= P.sum()
    te = 0.0
    Pb = P.sum(axis=0)                 # p(y0,x0)
    Pyb = P.sum(axis=2)               # p(y1,y0)
    Pbb = P.sum(axis=(0, 2))          # p(y0)
    for i in range(bins):
        for j in range(bins):
            for k in range(bins):
                p = P[i, j, k]
                if p <= 0:
                    continue
                p_y1_given_y0x0 = p / Pb[j, k]
                p_y1_given_y0 = Pyb[i, j] / Pbb[j] if Pbb[j] > 0 else 0
                if p_y1_given_y0 > 0:
                    te += p * np.log(p_y1_given_y0x0 / p_y1_given_y0)
    return float(te)


def partial_corr(a, b, z):
    """Partial correlation corr(a,b | z): remove the effect of a common driver z."""
    from numpy.linalg import lstsq
    z = np.asarray(z, float)
    Z = np.c_[np.ones_like(z), z] if z.ndim == 1 else np.c_[np.ones(len(z)), z]
    ra = a - Z.dot(lstsq(Z, a, rcond=None)[0])
    rb = b - Z.dot(lstsq(Z, b, rcond=None)[0])
    if ra.std() == 0 or rb.std() == 0:
        return 0.0
    return float(np.corrcoef(ra, rb)[0, 1])
