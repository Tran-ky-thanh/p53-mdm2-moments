# -*- coding: utf-8 -*-
"""
STEP 6 - Published Splatter noise + Pearson/MI/HSIC/dCor + permutation tests.
Results cached to data/cache_step6.npz; re-plotting only reloads.

Run:            python analysis/step6_splatter_noise_methods.py
Force recompute: RECOMPUTE=1 python analysis/step6_splatter_noise_methods.py
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
FIG = ROOT / "figures"; DATA = ROOT / "data"

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fast_ssa as F
import model as M
import association as A
import plotstyle
from cache import load_or_compute

N_CELLS = 3000
HORIZON = 900.0
DT = 15.0
T_SNAP = 600.0
DROPOUT_MIDS = [-1.0, 0.0, 1.0, 2.0, 3.0]
LIB_SCALE = 0.2
BCV = 0.5
N_REP = 5
PERM_MID = 1.0
CONDS = [(0.0, "CLOSED loop", 'C0'), (1.0, "NUTLIN-3", 'C3')]
METHODS = ['pearson', 'MI', 'HSIC', 'dCor']
MLAB = {'pearson': 'Pearson', 'MI': 'Mutual Info (nats)',
        'HSIC': 'normalized HSIC', 'dCor': 'distance corr'}


def compute():
    T = np.arange(0.0, HORIZON, DT)
    isnap = int(np.argmin(np.abs(T - T_SNAP)))
    snap = {}
    for nut, label, _ in CONDS:
        d = F.simulate_ensemble(N_CELLS, T, nutlin=nut, seed=6060 + int(nut * 10))
        snap[label] = np.stack([d[:, M.IDX['p53_mRNA'], isnap], d[:, M.IDX['Mdm2_mRNA'], isnap]], axis=1)
    out = {'mids': np.array(DROPOUT_MIDS), 'zerofrac': np.zeros(len(DROPOUT_MIDS))}
    for nut, label, _ in CONDS:
        C = snap[label]
        for m in METHODS:
            out[f'{label}__{m}__mean'] = np.zeros(len(DROPOUT_MIDS)); out[f'{label}__{m}__std'] = np.zeros(len(DROPOUT_MIDS))
        for di, mid in enumerate(DROPOUT_MIDS):
            vals = {m: [] for m in METHODS}; zf = []
            for r in range(N_REP):
                rng = np.random.default_rng(100 * di + r)
                obs = A.splatter_noise(C, rng, lib_scale=LIB_SCALE, bcv=BCV, dropout_mid=mid)
                zf.append(np.mean(obs == 0))
                meas = A.all_measures(obs[:, 0], obs[:, 1], n_kernel=1000, seed=r)
                for m in METHODS:
                    vals[m].append(meas[m])
            for m in METHODS:
                out[f'{label}__{m}__mean'][di] = np.nanmean(vals[m]); out[f'{label}__{m}__std'][di] = np.nanstd(vals[m])
            if nut == 0.0:
                out['zerofrac'][di] = np.mean(zf)
            print(f"[{label}] mid={mid:+.0f} | " + " ".join(f"{m}={out[f'{label}__{m}__mean'][di]:+.3f}" for m in METHODS))
        # permutation z at PERM_MID, averaged over several noise draws for stability
        zacc = {m: [] for m in METHODS}
        for rr in range(5):
            rng = np.random.default_rng(700 + rr)
            obs = A.splatter_noise(C, rng, lib_scale=LIB_SCALE, bcv=BCV, dropout_mid=PERM_MID)
            for m in METHODS:
                _, z, _ = A.permutation_z(obs[:, 0], obs[:, 1], method=m, B=200, n_sub=2000, seed=1)
                zacc[m].append(z)
        for m in METHODS:
            out[f'{label}__{m}__z'] = np.array([np.mean(zacc[m])])
            print(f"[{label}] perm {m}: z={np.mean(zacc[m]):+.1f} (+/-{np.std(zacc[m]):.1f})")
    return out


def plot(D):
    plotstyle.apply()
    xpos = np.arange(len(DROPOUT_MIDS))
    xlabels = [f"{m:+.0f}\n(~{z*100:.0f}% zero)" for m, z in zip(DROPOUT_MIDS, D['zerofrac'])]
    fig, axs = plt.subplots(1, 4, figsize=(19, 4.8))
    for ax, m in zip(axs, METHODS):
        for nut, label, col in CONDS:
            ax.errorbar(xpos, D[f'{label}__{m}__mean'], yerr=D[f'{label}__{m}__std'],
                        marker='o', color=col, capsize=3, label=label)
        ax.axhline(0, color='gray', lw=0.6)
        ax.set_xticks(xpos); ax.set_xticklabels(xlabels, fontsize=8)
        ax.set_title(MLAB[m], fontsize=11); ax.legend(fontsize=8); ax.grid(alpha=.3)
    fig.supxlabel("dropout_mid  (larger = more Splatter dropout)", fontsize=11, y=0.03)
    fig.suptitle(f"Step 6 - scRNA-seq technical noise (Splatter model, Zappia et al. 2017): detecting "
                 f"p53<->MDM2 mRNA dependence (Pearson/MI/HSIC/dCor)\n"
                 f"[N={N_CELLS} cells, snapshot t={T_SNAP:.0f} min, {N_REP} noise draws/point; "
                 f"HSIC/dCor subsampled to 1000]", fontsize=12)
    fig.tight_layout(rect=[0, 0.08, 1, 0.89])
    fig.savefig(str(FIG / "fig_step6_splatter_methods.png"), dpi=115)

    fig, ax = plt.subplots(figsize=(10, 5))
    w = 0.38; xpos = np.arange(len(METHODS))
    for k, (nut, label, col) in enumerate(CONDS):
        ax.bar(xpos + (k - 0.5) * w, [D[f'{label}__{m}__z'][0] for m in METHODS], width=w, color=col, label=label)
    ax.axhline(2, color='green', ls='--', lw=1, label='z=2 (~p<0.05)')
    ax.set_xticks(xpos); ax.set_xticklabels([MLAB[m] for m in METHODS], fontsize=9)
    ax.set_ylabel("permutation z-score (detection strength)")
    ax.set_title(f"Step 6 - Permutation z-scores at dropout_mid={PERM_MID:+.0f} "
                 f"(N={N_CELLS} cells, n_sub=2000, 200 permutations, mean over 5 noise draws)\n"
                 "CLOSED loop = truly dependent; NUTLIN = truly independent (z>0 there = library-size false positive)",
                 fontsize=10)
    ax.legend(); ax.grid(alpha=.3, axis='y')
    fig.tight_layout(); fig.savefig(str(FIG / "fig_step6_permutation_z.png"), dpi=120)
    print("Saved fig_step6_splatter_methods.png, fig_step6_permutation_z.png")


def main():
    D = load_or_compute(DATA / "cache_step6.npz", compute)
    plot(D)


if __name__ == "__main__":
    main()
