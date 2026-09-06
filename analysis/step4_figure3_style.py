# -*- coding: utf-8 -*-
"""
STEP 4 - Figure-3-style plots (Patterns 2021) for the pair (p53 mRNA, MDM2 mRNA).
Snapshot scatter + density contour + mean trajectory, plus moment time courses.
The original basal-state cache remains in data/cache_step4.npz.  This revised
Patterns-style run starts every dynamic species at zero and is saved separately
to data/cache_step4_zero_initial.npz; re-plotting only reloads that cache.

Run:            python analysis/step4_figure3_style.py
Force recompute: RECOMPUTE=1 python analysis/step4_figure3_style.py
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
FIG = ROOT / "figures"; DATA = ROOT / "data"

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
import fast_ssa as F
import model as M
import moments as mu
import plotstyle
from cache import load_or_compute

N_CELLS = 4000
HORIZON = 1200.0
DT = 4.0
T_SNAP = 300.0
CONDS = [(0.0, "CLOSED loop", 'C0'), (1.0, "NUTLIN-3", 'C3')]
KEYS = ['mean_p53', 'mean_mdm2', 'cov', 'corr', 'skew_p53', 'skew_mdm2']
CACHE = DATA / "cache_step4_zero_initial.npz"
INITIAL_STATE = np.zeros(len(M.species), dtype=float)


def compute():
    T = np.arange(0.0, HORIZON + 0.5 * DT, DT)
    isnap = int(np.argmin(np.abs(T - T_SNAP)))
    out = {'T': T, 'isnap': np.array([isnap]),
           'initial_state': INITIAL_STATE.copy(),
           'n_cells': np.array([N_CELLS]), 'dt': np.array([DT]),
           'horizon': np.array([HORIZON])}
    for nut, label, _ in CONDS:
        seed = 4040 + int(nut * 10)
        d = F.simulate_ensemble(N_CELLS, T, nutlin=nut,
                                x0=INITIAL_STATE, seed=seed)
        out[f'{label}__seed'] = np.array([seed])
        out[f'{label}__snapA'] = d[:, M.IDX['p53_mRNA'], isnap]
        out[f'{label}__snapB'] = d[:, M.IDX['Mdm2_mRNA'], isnap]
        s = mu.moment_time_courses(d[:, M.IDX['p53_mRNA'], :], d[:, M.IDX['Mdm2_mRNA'], :])
        for k in KEYS:
            out[f'{label}__{k}'] = s[k]
    return out


def smooth_series(values, window=11, nonnegative=False):
    """Savitzky-Golay display smoothing; cached statistics remain untouched."""
    values = np.asarray(values, float)
    out = values.copy()
    finite = np.flatnonzero(np.isfinite(values))
    if finite.size >= 5:
        width = min(window, finite.size if finite.size % 2 else finite.size - 1)
        width = max(5, width)
        out[finite] = savgol_filter(values[finite], width, 3, mode='interp')
        # Keep measured endpoints exact, especially the (0, 0) initial state.
        out[finite[0]] = values[finite[0]]
        out[finite[-1]] = values[finite[-1]]
    if nonnegative:
        out = np.maximum(out, 0.0)
    return out


def scatter_panel(ax, T, A_t, B_t, meanA, meanB, isnap, title, color):
    rng = np.random.default_rng(0)
    x = A_t + rng.uniform(-0.45, 0.45, A_t.shape); y = B_t + rng.uniform(-0.45, 0.45, B_t.shape)
    ax.scatter(x, y, s=4, alpha=0.10, color=color, edgecolors='none')
    xr = (max(0, A_t.min() - 2), A_t.max() + 2); yr = (max(0, B_t.min() - 2), B_t.max() + 2)
    Hh, xe, ye = np.histogram2d(A_t, B_t, bins=45, range=[xr, yr])
    Hh = gaussian_filter(Hh, sigma=1.15)
    positive = Hh[Hh > 0]
    levels = np.quantile(positive, [0.45, 0.60, 0.72, 0.82, 0.90])
    ax.contour(0.5*(xe[:-1]+xe[1:]), 0.5*(ye[:-1]+ye[1:]), Hh.T,
               levels=np.unique(levels), colors='k', linewidths=0.8, alpha=0.55)

    # Smooth only the displayed population-mean path. PCHIP draws a continuous
    # curve without overshooting between the denser 4-min raw snapshots.
    path_t = T[:isnap + 1]
    path_a = smooth_series(meanA[:isnap + 1], window=11, nonnegative=True)
    path_b = smooth_series(meanB[:isnap + 1], window=11, nonnegative=True)
    dense_t = np.linspace(path_t[0], path_t[-1], 500)
    dense_a = PchipInterpolator(path_t, path_a)(dense_t)
    dense_b = PchipInterpolator(path_t, path_b)(dense_t)
    ax.plot(dense_a, dense_b, '-', color='orange', lw=2.8, label='smoothed mean trajectory')
    ax.plot(meanA[0], meanB[0], 'o', color='green', ms=9, label='t=0')
    ax.plot(meanA[isnap], meanB[isnap], 's', color='red', ms=9, label=f't={T_SNAP:.0f}')
    ax.set_title(title, fontsize=13)
    ax.set_xlabel('p53 mRNA (counts)'); ax.set_ylabel('MDM2 mRNA (counts)'); ax.legend(fontsize=11)


def plot(D):
    plotstyle.apply()
    T = D['T']; isnap = int(D['isnap'][0])
    # 3 rows x 2 cols: row0 scatters, row1 mean+covariance, row2 skewness
    fig, axs = plt.subplots(3, 2, figsize=(13, 16))
    for j, (nut, label, col) in enumerate(CONDS):
        scatter_panel(axs[0, j], T, D[f'{label}__snapA'], D[f'{label}__snapB'],
                      D[f'{label}__mean_p53'], D[f'{label}__mean_mdm2'], isnap,
                      f"({'AB'[j]}) Snapshot t={T_SNAP:.0f} min - {label}\n"
                      f"cov={D[f'{label}__cov'][isnap]:+.1f}, corr={D[f'{label}__corr'][isnap]:+.2f}", col)
    axC = axs[1, 0]
    for nut, label, col in CONDS:
        lw = 3.0 if nut == 0.0 else 1.6
        axC.plot(T, smooth_series(D[f'{label}__mean_p53'], nonnegative=True),
                 color=col, ls='-', lw=lw, alpha=0.9, label=f'p53 mRNA - {label}')
        axC.plot(T, smooth_series(D[f'{label}__mean_mdm2'], nonnegative=True),
                 color=col, ls='--', lw=lw, alpha=0.9, label=f'MDM2 mRNA - {label}')
    axC.set_title("(C) Mean mRNA counts over time")
    axC.set_xlabel('time (min)'); axC.set_ylabel('mean count'); axC.legend(fontsize=10); axC.grid(alpha=.3)
    axC.text(0.5, 0.30, "solid p53-mRNA lines overlap and approach ~15:\np53 mRNA unchanged by Nutlin",
             transform=axC.transAxes, ha='center', fontsize=11, color='#555')
    labelmap = [(axs[1, 1], "(D) Covariance of p53 & MDM2 mRNA (across cells)", 'cov'),
                (axs[2, 0], "(E) Skewness of p53 mRNA (across cells)", 'skew_p53'),
                (axs[2, 1], "(F) Skewness of MDM2 mRNA (across cells)", 'skew_mdm2')]
    for ax, ti, key in labelmap:
        for nut, label, col in CONDS:
            ax.plot(T, smooth_series(D[f'{label}__{key}']),
                    color=col, label=label, lw=2.2)
        ax.axhline(0, color='gray', lw=0.6)
        ax.set_title(ti); ax.set_xlabel('time (min)')
        ax.set_ylabel(key.split('_')[0]); ax.legend(); ax.grid(alpha=.3)
    fig.suptitle(f"Step 4 - Figure-3-style plots (Patterns 2021) for the p53-MDM2 loop  [N={N_CELLS} cells]",
                 fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(str(FIG / "fig_step4_figure3_style.png"), dpi=120)
    print("Saved fig_step4_figure3_style.png")


def main():
    D = load_or_compute(CACHE, compute)
    plot(D)


if __name__ == "__main__":
    main()
