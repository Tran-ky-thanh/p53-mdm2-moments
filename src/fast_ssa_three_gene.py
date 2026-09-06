# -*- coding: utf-8 -*-
"""Numba Gillespie simulator for the six-species p53-MDM2-CDKN1A model."""
import numpy as np
from numba import njit, prange
import model_three_gene as M

PARAM_ORDER = [
    'k_txn_p53', 'k_deg_p53mRNA', 'k_tsl_p53', 'k_deg_p53', 'k_degp53',
    'Km_p53', 'k_txn_mdm2', 'Kd_mdm2', 'n_hill', 'k_txn_mdm2_0',
    'k_deg_mdm2mRNA', 'k_tsl_mdm2', 'k_deg_mdm2', 'k_txn_cdkn1a',
    'Kd_cdkn1a', 'n_cdkn1a', 'k_txn_cdkn1a_0', 'k_deg_cdkn1a_mRNA',
    'k_tsl_cdkn1a', 'k_deg_cdkn1a']

_STOICH = M.transitions.T.astype(np.int64)


def params_to_theta(params=None):
    p = dict(M.DEFAULT_PARAMS)
    if params:
        p.update(params)
    return np.array([p[k] for k in PARAM_ORDER], dtype=np.float64)


@njit(cache=True)
def _propensities(x, th, kdeg_eff):
    mp, P, md, D, mc, C = x
    hp = P**th[8] / (th[7]**th[8] + P**th[8])
    hc = P**th[15] / (th[14]**th[15] + P**th[15])
    a = np.empty(13)
    a[0] = th[0]; a[1] = th[1]*mp; a[2] = th[2]*mp; a[3] = th[3]*P
    a[4] = kdeg_eff*D*P/(th[5]+P)
    a[5] = th[6]*hp+th[9]; a[6] = th[10]*md; a[7] = th[11]*md; a[8] = th[12]*D
    a[9] = th[13]*hc+th[16]; a[10] = th[17]*mc; a[11] = th[18]*mc; a[12] = th[19]*C
    return a


@njit(cache=True)
def _simulate_cell(th, nutlin, x0, times, seed):
    np.random.seed(seed)
    out = np.empty((x0.size, times.size), dtype=np.float64)
    x = x0.astype(np.float64).copy(); out[:, 0] = x; t = times[0]
    kdeg_eff = th[4] * (1.0 - nutlin)
    for i in range(1, times.size):
        target = times[i]
        while True:
            a = _propensities(x, th, kdeg_eff); a0 = a.sum()
            if a0 <= 0.0:
                t = target; break
            r1 = np.random.random()
            while r1 == 0.0:
                r1 = np.random.random()
            tau = np.log(1.0/r1)/a0
            if t + tau > target:
                t = target; break
            t += tau; r2 = np.random.random()*a0; acc = 0.0; chosen = 12
            for j in range(13):
                acc += a[j]
                if acc > r2:
                    chosen = j; break
            x += _STOICH[chosen]
        out[:, i] = x
    return out


@njit(parallel=True, cache=True)
def _simulate_ensemble(th, nutlin, x0, times, n, seed):
    out = np.empty((n, x0.size, times.size), dtype=np.float64)
    for cell in prange(n):
        out[cell] = _simulate_cell(th, nutlin, x0, times, seed + cell)
    return out


def simulate_ensemble(n, times, params=None, nutlin=0.0, x0=None, seed=0):
    """Simulate independent cells and retain only the requested observation times."""
    if x0 is None:
        x0 = M.initial_state
    return _simulate_ensemble(params_to_theta(params), float(nutlin),
                              np.asarray(x0, float), np.asarray(times, float), int(n),
                              int(seed) * 1000003)
