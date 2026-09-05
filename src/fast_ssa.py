# -*- coding: utf-8 -*-
"""
numba-accelerated (JIT) Gillespie simulator, SPECIALIZED to the p53-MDM2 model (MM/Hill)
defined in model.py. Used to simulate LARGE cell POPULATIONS (scRNA-seq).

IMPORTANT: the propensities and stoichiometry here must MATCH model.py EXACTLY (the reference
that runs on the Patterns 2021 ssa.py engine). step1b_validate_fast_ssa.py checks that both give
identical propensities/statistics -> equivalent.

Order of the parameters in the `theta` array (matching model.DEFAULT_PARAMS):
 0:k_txn_p53 1:k_deg_p53mRNA 2:k_tsl_p53 3:k_deg_p53 4:k_degp53 5:Km_p53
 6:k_txn_mdm2 7:Kd_mdm2 8:n_hill 9:k_txn_mdm2_0 10:k_deg_mdm2mRNA 11:k_tsl_mdm2 12:k_deg_mdm2
State x = [p53_mRNA, p53, Mdm2_mRNA, Mdm2].
"""

import numpy as np
from numba import njit, prange
import model as M

PARAM_ORDER = ['k_txn_p53', 'k_deg_p53mRNA', 'k_tsl_p53', 'k_deg_p53', 'k_degp53',
               'Km_p53', 'k_txn_mdm2', 'Kd_mdm2', 'n_hill', 'k_txn_mdm2_0',
               'k_deg_mdm2mRNA', 'k_tsl_mdm2', 'k_deg_mdm2']


def params_to_theta(params=None):
    p = dict(M.DEFAULT_PARAMS)
    if params:
        p.update(params)
    return np.array([p[k] for k in PARAM_ORDER], dtype=np.float64)


@njit(cache=True)
def _propensities(x, th, kdeg_eff):
    mp = x[0]; P = x[1]; md = x[2]; D = x[3]
    hill = P**th[8] / (th[7]**th[8] + P**th[8])
    mm = D * P / (th[5] + P)
    a = np.empty(9)
    a[0] = th[0]                       # R0 p53 transcription (constant)
    a[1] = th[1] * mp                  # R1 p53 mRNA degradation
    a[2] = th[2] * mp                  # R2 p53 translation
    a[3] = th[3] * P                   # R3 basal p53 degradation
    a[4] = kdeg_eff * mm               # R4 MDM2-mediated p53 degradation (MM) [Nutlin]
    a[5] = th[6] * hill + th[9]        # R5 p53-activated MDM2 transcription (Hill)
    a[6] = th[10] * md                 # R6 MDM2 mRNA degradation
    a[7] = th[11] * md                 # R7 MDM2 translation
    a[8] = th[12] * D                  # R8 MDM2 degradation
    return a


# stoichiometry (9 reactions x 4 species), matching model.transitions
_STOICH = np.array([
    ( 1, 0, 0, 0),  # R0
    (-1, 0, 0, 0),  # R1
    ( 0, 1, 0, 0),  # R2
    ( 0,-1, 0, 0),  # R3
    ( 0,-1, 0, 0),  # R4
    ( 0, 0, 1, 0),  # R5
    ( 0, 0,-1, 0),  # R6
    ( 0, 0, 0, 1),  # R7
    ( 0, 0, 0,-1),  # R8
], dtype=np.int64)


@njit(cache=True)
def _simulate_cell(th, nutlin, x0, T, stoich, seed):
    np.random.seed(seed)               # per-cell seed -> reproducible regardless of thread
    n_sp = x0.shape[0]
    n_t = T.shape[0]
    out = np.empty((n_sp, n_t), dtype=np.float64)
    x = x0.astype(np.float64).copy()
    kdeg_eff = th[4] * (1.0 - nutlin)
    out[:, 0] = x
    t = T[0]
    for i in range(1, n_t):
        t_target = T[i]
        while True:
            a = _propensities(x, th, kdeg_eff)
            a0 = 0.0
            for j in range(9):
                a0 += a[j]
            if a0 <= 0.0:
                t = t_target
                break
            r1 = np.random.random()
            while r1 == 0.0:
                r1 = np.random.random()
            tau = (1.0 / a0) * np.log(1.0 / r1)
            if t + tau > t_target:
                t = t_target
                break
            t += tau
            # choose reaction
            r2 = np.random.random() * a0
            cum = 0.0
            jsel = 8
            for j in range(9):
                cum += a[j]
                if cum > r2:
                    jsel = j
                    break
            for k in range(n_sp):
                x[k] += stoich[jsel, k]
        out[:, i] = x
    return out


@njit(parallel=True, cache=True)
def _simulate_ensemble(th, nutlin, x0, T, N, stoich, base_seed):
    n_sp = x0.shape[0]
    n_t = T.shape[0]
    data = np.empty((N, n_sp, n_t), dtype=np.float64)
    for c in prange(N):
        data[c] = _simulate_cell(th, nutlin, x0, T, stoich, base_seed + c)
    return data


def simulate_ensemble(N, T, params=None, nutlin=0.0, x0=None, seed=None):
    """
    Simulate N independent cells; return an array (N, n_species, n_time).
    Equivalent to ssa.SSA_Fixed_Width_Trajectory but JIT-compiled and parallel.
    """
    base_seed = 0 if seed is None else int(seed)
    th = params_to_theta(params)
    if x0 is None:
        x0 = M.initial_state.copy()
    x0 = np.asarray(x0, dtype=np.float64)
    T = np.asarray(T, dtype=np.float64)
    # per-cell seeds base_seed*N_off + c make the whole ensemble reproducible across runs/threads
    return _simulate_ensemble(th, float(nutlin), x0, T, int(N), _STOICH, base_seed * 1000003)
