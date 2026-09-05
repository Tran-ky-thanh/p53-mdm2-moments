# -*- coding: utf-8 -*-
"""
STEP 3 - Higher-order MOMENTS from "simulated scRNA-seq" (core Patterns 2021 idea).
Mean/variance/covariance/skewness of the (p53 mRNA, MDM2 mRNA) pair, closed vs Nutlin.
Simulation cached to data/cache_step3.npz; re-plotting only reloads.

Run:            python analysis/step3_scRNAseq_moments.py
Force recompute: RECOMPUTE=1 python analysis/step3_scRNAseq_moments.py
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
import moments as mu
import plotstyle
from cache import load_or_compute

N_CELLS = 3000
HORIZON = 1500.0
DT = 15.0
CONDS = [(0.0, "CLOSED loop (no drug)", 'C0'), (1.0, "NUTLIN-3 (open loop)", 'C3')]
KEYS = ['mean_p53', 'mean_mdm2', 'var_p53', 'var_mdm2', 'cov', 'corr', 'skew_p53', 'skew_mdm2']


def compute():
    T = np.arange(0.0, HORIZON, DT)
    out = {'T': T}
    for nut, label, _ in CONDS:
        data = F.simulate_ensemble(N_CELLS, T, nutlin=nut, seed=2024 + int(nut * 10))
        s = mu.moment_time_courses(data[:, M.IDX['p53_mRNA'], :], data[:, M.IDX['Mdm2_mRNA'], :])
        for k in KEYS:
            out[f'{label}__{k}'] = s[k]
        print(f"[{label}] late cov={s['cov'][-1]:+.1f} skew(MDM2mRNA)={s['skew_mdm2'][-1]:+.2f}")
    return out


def plot(D):
    plotstyle.apply()
    T = D['T']
    # 3 rows x 2 cols so each panel is large and readable
    panels = [('mean_p53', "Mean of p53 mRNA (per cell)", 'count'),
              ('mean_mdm2', "Mean of MDM2 mRNA (per cell)", 'count'),
              ('cov', "Covariance of p53 & MDM2 mRNA (across cells)", 'covariance'),
              ('var_p53', "Variance of p53 mRNA (across cells)", 'variance'),
              ('skew_p53', "Skewness of p53 mRNA (across cells)", 'skewness'),
              ('skew_mdm2', "Skewness of MDM2 mRNA (across cells)", 'skewness')]
    fig, axs = plt.subplots(3, 2, figsize=(13, 15))
    for ax, (key, ti, yl) in zip(axs.ravel(), panels):
        for nut, label, col in CONDS:
            lw = 3.0 if nut == 0.0 else 1.6     # CLOSED thick under, NUTLIN thin on top
            ax.plot(T, D[f'{label}__{key}'], color=col, label=label, lw=lw, alpha=0.92)
        ax.set_title(ti); ax.set_ylabel(yl); ax.set_xlabel('time (min)')
        ax.grid(alpha=0.3); ax.legend()
    axs[0, 0].text(0.5, 0.08, "the two lines overlap:\np53 mRNA is the same in both conditions",
                   transform=axs[0, 0].transAxes, ha='center', fontsize=11, color='#555')
    fig.suptitle(f"Step 3 - Population statistics of the (p53, MDM2) mRNA pair over time  "
                 f"[N={N_CELLS} cells]\nmean (1st moment), covariance (2nd), skewness (3rd)", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(str(FIG / "fig_step3_moments.png"), dpi=120)
    print("Saved fig_step3_moments.png")


def main():
    D = load_or_compute(DATA / "cache_step3.npz", compute)
    plot(D)


if __name__ == "__main__":
    main()
