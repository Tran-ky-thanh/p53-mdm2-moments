# -*- coding: utf-8 -*-
"""
STEP 5 - Simple dropout noise + Pearson/MI/HSIC/dCor between p53 and MDM2 mRNA.
Results cached to data/cache_step5.npz; re-plotting only reloads.

Run:            python analysis/step5_noise_association_methods.py
Force recompute: RECOMPUTE=1 python analysis/step5_noise_association_methods.py
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

N_CELLS = 4000
HORIZON = 900.0
DT = 15.0
T_SNAP = 600.0
BETAS = [1.0, 0.6, 0.3, 0.15, 0.08]
N_REP = 6
N_KERNEL = 1200
CONDS = [(0.0, "CLOSED loop", 'C0'), (1.0, "NUTLIN-3", 'C3')]
METHODS = ['pearson', 'MI', 'HSIC', 'dCor']
METHOD_LABEL = {'pearson': 'Pearson correlation',
                'MI': 'Mutual information (nats)', 'HSIC': 'normalized HSIC',
                'dCor': 'distance correlation'}


def compute():
    T = np.arange(0.0, HORIZON, DT)
    isnap = int(np.argmin(np.abs(T - T_SNAP)))
    out = {'betas': np.array(BETAS), 'n_cells': np.array([N_CELLS]),
           't_snap': np.array([T_SNAP]), 'n_rep': np.array([N_REP]),
           'n_kernel': np.array([N_KERNEL])}
    for nut, label, _ in CONDS:
        d = F.simulate_ensemble(N_CELLS, T, nutlin=nut, seed=5050 + int(nut * 10))
        Xc = d[:, M.IDX['p53_mRNA'], isnap]; Yc = d[:, M.IDX['Mdm2_mRNA'], isnap]
        # Retain the clean snapshot so future metric/plot additions do not need
        # to rerun the Gillespie trajectories.
        out[f'{label}__clean_p53_mRNA'] = Xc
        out[f'{label}__clean_mdm2_mRNA'] = Yc
        for m in METHODS:
            out[f'{label}__{m}__mean'] = np.zeros(len(BETAS)); out[f'{label}__{m}__std'] = np.zeros(len(BETAS))
        for bi, beta in enumerate(BETAS):
            vals = {m: [] for m in METHODS}
            for r in range(N_REP):
                rng = np.random.default_rng(1000 * bi + r)
                meas = A.all_measures(A.add_dropout(Xc, beta, rng), A.add_dropout(Yc, beta, rng),
                                      n_kernel=N_KERNEL, seed=r)
                for m in METHODS:
                    vals[m].append(meas[m])
            for m in METHODS:
                out[f'{label}__{m}__mean'][bi] = np.mean(vals[m]); out[f'{label}__{m}__std'][bi] = np.std(vals[m])
            print(f"[{label}] beta={beta:.2f} | " + " ".join(f"{m}={out[f'{label}__{m}__mean'][bi]:+.3f}" for m in METHODS))
    return out


def plot(D):
    plotstyle.apply()
    fig, axes = plt.subplots(1, 4, figsize=(19, 5.0))
    xpos = np.arange(len(BETAS))
    for ax, m in zip(axes, METHODS):
        for nut, label, col in CONDS:
            ax.errorbar(xpos, D[f'{label}__{m}__mean'], yerr=D[f'{label}__{m}__std'],
                        marker='o', lw=2.3, ms=6, color=col, capsize=3, label=label)
        ax.axhline(0, color='gray', lw=0.6)
        ax.set_xticks(xpos); ax.set_xticklabels([f"{b:.2f}" for b in BETAS])
        ax.set_title(METHOD_LABEL[m], fontsize=13); ax.legend(fontsize=10); ax.grid(alpha=.3)
    fig.supxlabel("capture efficiency beta  (lower = more dropout)", fontsize=12, y=0.03)
    fig.suptitle(f"Step 5 - Detecting p53<->MDM2 mRNA dependence under simple dropout noise  "
                 f"[N={N_CELLS} cells, snapshot t={T_SNAP:.0f} min, {N_REP} noise draws/point]\n"
                 "(closed loop = dependent; Nutlin = open loop, dependence ~0)", fontsize=12)
    fig.tight_layout(rect=[0, 0.08, 1, 0.89])
    fig.savefig(str(FIG / "fig_step5_noise_association.png"), dpi=120)
    print("Saved fig_step5_noise_association.png")


def main():
    cache_path = DATA / "cache_step5.npz"
    D = load_or_compute(cache_path, compute)
    metadata = {'n_cells': np.array([N_CELLS]), 't_snap': np.array([T_SNAP]),
                'n_rep': np.array([N_REP]), 'n_kernel': np.array([N_KERNEL])}
    if any(key not in D for key in metadata):
        D.update(metadata)
        np.savez_compressed(cache_path, **D)
        print(f"[cache] added self-describing metadata to {cache_path.name}")
    plot(D)


if __name__ == "__main__":
    main()
