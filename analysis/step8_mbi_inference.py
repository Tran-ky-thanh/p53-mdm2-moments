# -*- coding: utf-8 -*-
"""
STEP 8 - Run THEIR method: NON-LINEAR MBI (moment-based inference), vendored in src/mbi,
on the simulated moment time-courses of the mRNA pair p53-MDM2.
Recovers the directed edge p53->MDM2 (A->B) in the closed loop.
The 8000-cell example fit is cached to data/cache_step8.npz. A separate
multi-seed robustness analysis (raw moments, fitted rates and networks) is
cached to data/cache_step8_robust.npz; re-plotting reloads both caches.

Run:            python analysis/step8_mbi_inference.py
Force recompute: RECOMPUTE=1 python analysis/step8_mbi_inference.py
"""
import sys, pathlib
import contextlib
import io
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
FIG = ROOT / "figures"; DATA = ROOT / "data"

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import nnls
import fast_ssa as F
import model as M
import plotstyle
from cache import load_or_compute
from mbi.two_Dim_mRNA import generate_design_blocks as GDB
from mbi.NLLS import NLLS_Fit
from mbi.NLLS_Push_Forward import PushForward_Func
from mbi.Regulatory_Rules import Give_Regulatory_Network

hash_lookUp = {(0, 0): 0, (1, 0): 1, (0, 1): 2, (2, 0): 3, (1, 1): 4, (0, 2): 5,
               (3, 0): 6, (2, 1): 7, (1, 2): 8, (0, 3): 9,
               (4, 0): 10, (3, 1): 11, (2, 2): 12, (1, 3): 13, (0, 4): 14}
moms_lookUp_all = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2),
                   (3, 0), (2, 1), (1, 2), (0, 3), (4, 0), (3, 1), (2, 2), (1, 3), (0, 4)]
NLLS_Moments_Fit = np.array([True]*6 + [False]*4 + [False]*5)
Moments_Spline = np.array([False]*6 + [True]*4 + [False]*5)
Spline_der_bool = True
Mom_inds = [0, 0, 1, 2, 4, 4, 2, 1, 2, 1]
N_CELLS = 8000
HORIZON = 1000.0
DT = 2.0
SHIFT = 10
SUBSAMPLE = 10
CONDS = [(0.0, "CLOSED loop", 'C0'), (1.0, "NUTLIN-3", 'C3')]
ROBUST_N_CELLS = 2500
ROBUST_SEEDS = np.array([11, 23, 37, 51, 67], dtype=int)
ROBUST_SUBSAMPLE = 20


def raw_moments(dataA, dataB):
    n_t = dataA.shape[1]; E = np.zeros((n_t, 15))
    for (l1, l2), idx in hash_lookUp.items():
        E[:, idx] = np.mean((dataA ** l1) * (dataB ** l2), axis=0)
    return E


def linear_init(E, T, DB):
    fit_rows = np.where(NLLS_Moments_Fit[:DB.shape[1]])[0]
    dEdt = np.gradient(E[:, fit_rows], T, axis=0)
    n_t = len(T); nK = DB.shape[0]
    R = np.zeros((n_t * len(fit_rows), nK)); b = np.zeros(n_t * len(fit_rows)); r = 0
    for ti in range(n_t):
        for j, row in enumerate(fit_rows):
            R[r, :] = DB[:, row, :].dot(E[ti, :]); b[r] = dEdt[ti, j]; r += 1
    return nnls(R, b)[0]


def infer_one(nut, seed, n_cells, subsample, quiet=False):
    DB = GDB(hash_lookUp, moms_lookUp_all[:10], Verbose=False)
    W = np.eye(int(NLLS_Moments_Fit.sum())); W[1, 1] = 40.0; W[2, 2] = 40.0
    T_full = np.arange(0.0, HORIZON, DT)
    d = F.simulate_ensemble(n_cells, T_full, nutlin=nut, seed=int(seed))
    E = raw_moments(d[:, M.IDX['p53_mRNA'], :], d[:, M.IDX['Mdm2_mRNA'], :])[SHIFT::subsample, :]
    T = T_full[SHIFT::subsample]
    init = np.abs(linear_init(E, T, DB))
    if quiet:
        with contextlib.redirect_stdout(io.StringIO()):
            K, _ = NLLS_Fit(init, E, T, DB, NLLS_Moments_Fit, Moments_Spline,
                            Spline_der_bool, weights=W)
    else:
        K, _ = NLLS_Fit(init, E, T, DB, NLLS_Moments_Fit, Moments_Spline,
                        Spline_der_bool, weights=W)
    react = np.array([E[:, Mom_inds[k]] * K[k] for k in range(len(K))])
    reg = Give_Regulatory_Network(T[-1] - T[0], np.average(react, axis=1))
    _, pred = PushForward_Func(K, E, T, DB, NLLS_Moments_Fit, Moments_Spline,
                              Spline_der_bool)
    return T, E, pred, K, reg


def compute():
    out = {}
    for nut, label, _ in CONDS:
        T, E, pred, K, reg = infer_one(nut, 808 + int(nut * 10), N_CELLS, SUBSAMPLE)
        out[f'{label}__T'] = T; out[f'{label}__E'] = E; out[f'{label}__pred'] = pred
        out[f'{label}__K'] = K; out[f'{label}__reg'] = reg
        print(f"[{label}] A->B (p53->MDM2)={reg[1,0]:+.2f}  B->A={reg[0,1]:+.2f}")
    return out


def compute_robust():
    out = {'seeds': ROBUST_SEEDS, 'n_cells': np.array([ROBUST_N_CELLS]),
           'subsample': np.array([ROBUST_SUBSAMPLE])}
    for nut, label, _ in CONDS:
        moments, rates, networks = [], [], []
        for seed in ROBUST_SEEDS:
            T, E, _, K, reg = infer_one(nut, seed, ROBUST_N_CELLS,
                                        ROBUST_SUBSAMPLE, quiet=True)
            moments.append(E); rates.append(K); networks.append(reg)
            print(f"[robust {label}] seed={seed}: p53->MDM2={reg[1,0]:+.3f}, "
                  f"MDM2->p53={reg[0,1]:+.3f}")
        out[f'{label}__T'] = T
        out[f'{label}__E_reps'] = np.asarray(moments)
        out[f'{label}__K_reps'] = np.asarray(rates)
        out[f'{label}__reg_reps'] = np.asarray(networks)
    return out


def plot(D, R):
    plotstyle.apply()
    fig, axs = plt.subplots(3, 2, figsize=(13, 16))
    for col, (nut, label, cc) in enumerate(CONDS):
        T = D[f'{label}__T']; E = D[f'{label}__E']; pred = D[f'{label}__pred']; reg = D[f'{label}__reg']
        ax = axs[0, col]
        ax.plot(T, E[:, 1], 'o', ms=4, color='C0', label='E[A] data'); ax.plot(T, pred[:, 1], '-', lw=2, color='C0', label='E[A] MBI fit')
        ax.plot(T, E[:, 2], 'o', ms=4, color='C1', label='E[B] data'); ax.plot(T, pred[:, 2], '-', lw=2, color='C1', label='E[B] MBI fit')
        ax.set_title(f"({'A' if col == 0 else 'B'}) {label}: mean fit (NLLS push-forward)")
        ax.set_xlabel('time (min)'); ax.set_ylabel('moment'); ax.legend(fontsize=11)

        reps = R[f'{label}__reg_reps']
        reg_mean = reps.mean(axis=0); reg_sd = reps.std(axis=0, ddof=1)
        ax = axs[1, col]
        all_means = [R[f'{lab}__reg_reps'].mean(axis=0) for _, lab, _ in CONDS]
        vlim = max(1, max(np.abs(m).max() for m in all_means))
        im = ax.imshow(reg_mean, cmap='RdBu_r', vmin=-vlim, vmax=vlim)
        ax.set_xticks([0, 1]); ax.set_xticklabels(['p53(A)', 'MDM2(B)'])
        ax.set_yticks([0, 1]); ax.set_yticklabels(['p53(A)', 'MDM2(B)'])
        ax.set_xlabel('REGULATOR'); ax.set_ylabel('REGULATED')
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{reg_mean[i,j]:+.2f}\n+/- {reg_sd[i,j]:.2f}",
                        ha='center', va='center', fontsize=13, fontweight='bold')
        panel = 'C' if col == 0 else 'D'
        ax.set_title(f"({panel}) {label}: mean MBI network +/- SD\n"
                     f"({len(ROBUST_SEEDS)} seeds x {ROBUST_N_CELLS} cells)")
        fig.colorbar(im, ax=ax, fraction=0.046)
    axc = axs[2, 0]; labels = [c[1] for c in CONDS]
    ab_reps = [R[f'{l}__reg_reps'][:, 1, 0] for l in labels]
    ba_reps = [R[f'{l}__reg_reps'][:, 0, 1] for l in labels]
    ab = [x.mean() for x in ab_reps]; ba = [x.mean() for x in ba_reps]
    ab_sd = [x.std(ddof=1) for x in ab_reps]; ba_sd = [x.std(ddof=1) for x in ba_reps]
    x = np.arange(2); w = 0.35
    axc.bar(x - w/2, ab, w, yerr=ab_sd, capsize=5,
            label='A->B  (p53->MDM2)', color='C2')
    axc.bar(x + w/2, ba, w, yerr=ba_sd, capsize=5,
            label='B->A  (MDM2->p53)', color='C4')
    axc.axhline(0, color='k', lw=0.6); axc.set_xticks(x); axc.set_xticklabels(labels)
    axc.set_ylabel('regulation strength (flux)'); axc.legend(fontsize=11)
    axc.set_title(f"(E) Inferred edges: mean +/- SD across {len(ROBUST_SEEDS)} seeds")
    axs[2, 1].axis('off')
    fig.suptitle(f"Step 8 - Non-linear moment-based inference (Raharinirina et al. 2021) on simulated p53-MDM2 moments\n"
                 f"[A-B: N={N_CELLS} cells/condition; C-E: {len(ROBUST_SEEDS)} seeds x "
                 f"{ROBUST_N_CELLS} cells; moments up to order 3-4]", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(str(FIG / "fig_step8_mbi.png"), dpi=120)
    print("Saved fig_step8_mbi.png")


def main():
    D = load_or_compute(DATA / "cache_step8.npz", compute)
    R = load_or_compute(DATA / "cache_step8_robust.npz", compute_robust)
    plot(D, R)


if __name__ == "__main__":
    main()
