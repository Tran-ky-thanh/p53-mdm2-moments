# -*- coding: utf-8 -*-
"""
STEP 6b - Does standard scRNA-seq PREPROCESSING (QC filtering, library-size
normalization, log transform, scaling) remove the spurious dependence that
Step 6 found under Splatter noise?

Step 6 fed RAW noisy counts straight into Pearson/MI/HSIC/dCor. Real pipelines
(Seurat, Scanpy) never do that. Here we re-run the same snapshot through five
preprocessing pipelines and compare each measure against the CLEAN-count ground
truth.

Key design point: library-size normalization is only meaningful with a
transcriptome-wide total, so alongside the two genes of interest we simulate a
background of N_BG independent "spectator" genes that share the same per-cell
library-size factor - exactly the confounder present in real data.

Run:            python analysis/step6b_normalization.py
Force recompute: RECOMPUTE=1 python analysis/step6b_normalization.py
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

# --- kept identical to Step 6 so the comparison is apples-to-apples ---
N_CELLS = 3000
HORIZON = 900.0
DT = 15.0
T_SNAP = 600.0
LIB_SCALE = 0.2
BCV = 0.5
DROPOUT_MID = 1.0
CONDS = [(0.0, "CLOSED loop", 'C0'), (1.0, "NUTLIN-3", 'C3')]
METHODS = ['pearson', 'MI', 'HSIC', 'dCor']
MLAB = {'pearson': 'Pearson r', 'MI': 'Mutual Info (nats)',
        'HSIC': 'normalized HSIC', 'dCor': 'distance correlation'}

# --- background transcriptome (needed for a real library size) ---
N_BG = 300
BG_MEAN_LOG = np.log(5.0)
BG_MEAN_SIG = 1.0

N_REP = 5          # noise draws for the point estimates
N_REP_Z = 3        # noise draws for the permutation z-scores
PERM_B = 150
PERM_NSUB = 1200

PIPELINES = [
    ('raw',      'raw counts\n(no preproc.)'),
    ('log',      'log1p\nonly'),
    ('cp10k',    'CP10k\n+ log1p'),
    ('cp10k_qc', 'QC +\nCP10k + log1p'),
    ('resid',    'CP10k + log1p\n+ regress depth'),
]


def make_background(n_cells, rng):
    """Independent spectator genes: per-gene mean is fixed, counts are Poisson."""
    mu = rng.lognormal(mean=BG_MEAN_LOG, sigma=BG_MEAN_SIG, size=(1, N_BG))
    return rng.poisson(np.repeat(mu, n_cells, axis=0)).astype(float)


def preprocess(obs, kind):
    """
    obs : (n_cells, 2 + N_BG) observed counts, columns 0/1 = p53 and MDM2 mRNA.
    Returns (x, y) ready for the dependence measures.
    """
    tot = obs.sum(axis=1)
    keep = np.ones(obs.shape[0], bool)

    if kind == 'cp10k_qc':
        # typical QC: drop the lowest-depth cells and cells with too few genes detected
        n_det = (obs > 0).sum(axis=1)
        keep = (tot >= np.percentile(tot, 5)) & (n_det >= np.percentile(n_det, 5))
        obs = obs[keep]; tot = tot[keep]

    g = obs[:, :2].astype(float)

    if kind == 'raw':
        return g[:, 0], g[:, 1]
    if kind == 'log':
        z = np.log1p(g)
        return z[:, 0], z[:, 1]

    # library-size normalization: counts per 10,000, then log1p
    f = 1e4 / np.maximum(tot, 1.0)
    z = np.log1p(g * f[:, None])
    if kind in ('cp10k', 'cp10k_qc'):
        return z[:, 0], z[:, 1]

    if kind == 'resid':
        # residualize each gene on log sequencing depth (removes any leftover
        # depth-driven component that CP10k did not fully absorb)
        d = np.log(np.maximum(tot, 1.0))
        Xd = np.stack([np.ones_like(d), d], axis=1)
        beta, *_ = np.linalg.lstsq(Xd, z, rcond=None)
        r = z - Xd @ beta
        return r[:, 0], r[:, 1]

    raise ValueError(kind)


N_SWEEP = [1500, 3000, 6000, 12000]
N_MAX = max(N_SWEEP)
N_REP_SWEEP = 3


def compute_ncells():
    """
    How does significance scale with the number of cells?

    Uses the Fisher transform, which gives the significance of a Pearson
    correlation analytically (no permutations needed):
        z = arctanh(r) * sqrt(n - 3).
    One N_MAX-cell ensemble is simulated per condition and subsampled, so the
    curves differ only in sample size, not in the underlying population.
    """
    T = np.arange(0.0, HORIZON, DT)
    isnap = int(np.argmin(np.abs(T - T_SNAP)))
    out = {'n_sweep': np.array(N_SWEEP, float)}
    for nut, label, _ in CONDS:
        d = F.simulate_ensemble(N_MAX, T, nutlin=nut, seed=6060 + int(nut * 10))
        C = np.stack([d[:, M.IDX['p53_mRNA'], isnap],
                      d[:, M.IDX['Mdm2_mRNA'], isnap]], axis=1)
        for p in ('raw', 'cp10k'):
            out[f'{label}__{p}__r'] = np.zeros(len(N_SWEEP))
            out[f'{label}__{p}__zf'] = np.zeros(len(N_SWEEP))
        for ni, n in enumerate(N_SWEEP):
            acc = {p: [] for p in ('raw', 'cp10k')}
            for r in range(N_REP_SWEEP):
                rng = np.random.default_rng(5500 + 23 * r)
                sub = rng.choice(N_MAX, n, replace=False)
                full = np.concatenate([C[sub], make_background(n, rng)], axis=1)
                obs = A.splatter_noise(full, rng, lib_scale=LIB_SCALE, bcv=BCV,
                                       dropout_mid=DROPOUT_MID)
                for p in ('raw', 'cp10k'):
                    x, y = preprocess(obs, p)
                    acc[p].append(A.pearson(x, y))
            for p in ('raw', 'cp10k'):
                rbar = float(np.mean(acc[p]))
                out[f'{label}__{p}__r'][ni] = rbar
                out[f'{label}__{p}__zf'][ni] = np.arctanh(rbar) * np.sqrt(n - 3)
            print(f"[{label}] n={n:6d} | raw r={out[f'{label}__raw__r'][ni]:+.4f} "
                  f"z={out[f'{label}__raw__zf'][ni]:+6.1f} | "
                  f"cp10k r={out[f'{label}__cp10k__r'][ni]:+.4f} "
                  f"z={out[f'{label}__cp10k__zf'][ni]:+6.1f}")
    return out


def plot_ncells(D):
    plotstyle.apply()
    n = D['n_sweep']
    fig, ax = plt.subplots(figsize=(10, 7))
    style = {'raw': ('-', 'o', 'raw counts'), 'cp10k': ('--', 's', 'CP10k + log1p')}
    for _, label, col in CONDS:
        for p, (ls, mk, plab) in style.items():
            ax.plot(n, D[f'{label}__{p}__zf'], ls, marker=mk, ms=9, lw=2.2, color=col,
                    fillstyle='full' if p == 'raw' else 'none',
                    label=f"{label} - {plab}")
    ax.axhline(2, color='green', ls=':', lw=1.5, label='z = 2  (~p < 0.05)')
    ax.axhline(0, color='k', lw=0.6)
    ax.set_xscale('log'); ax.set_xticks(N_SWEEP)
    ax.set_xticklabels([str(int(v)) for v in N_SWEEP])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("number of cells")
    ax.set_ylabel("Fisher z of Pearson r   (significance)")
    ax.set_title("More cells make the library-size artefact MORE significant\n"
                 "NUTLIN-3 is truly independent: its raw-count curve is pure false positive",
                 fontsize=14)
    ax.legend(fontsize=11); ax.grid(alpha=.3, which='both')
    fig.tight_layout()
    fig.savefig(str(FIG / "fig_step6b_ncells.png"), dpi=120)
    print("Saved fig_step6b_ncells.png")


def compute():
    T = np.arange(0.0, HORIZON, DT)
    isnap = int(np.argmin(np.abs(T - T_SNAP)))
    out = {'n_cells': np.array([N_CELLS]), 'n_bg': np.array([N_BG])}

    for nut, label, _ in CONDS:
        d = F.simulate_ensemble(N_CELLS, T, nutlin=nut, seed=6060 + int(nut * 10))
        C = np.stack([d[:, M.IDX['p53_mRNA'], isnap],
                      d[:, M.IDX['Mdm2_mRNA'], isnap]], axis=1)

        # ---- ground truth on CLEAN counts (no technical noise at all) ----
        clean = A.all_measures(C[:, 0], C[:, 1], n_kernel=1000, seed=0)
        for m in METHODS:
            out[f'{label}__clean__{m}'] = np.array([clean[m]])
        print(f"[{label}] CLEAN truth: " + " ".join(f"{m}={clean[m]:+.3f}" for m in METHODS))

        acc = {(p, m): [] for p, _ in PIPELINES for m in METHODS}
        for r in range(N_REP):
            rng = np.random.default_rng(4200 + 17 * r)
            full = np.concatenate([C, make_background(N_CELLS, rng)], axis=1)
            obs = A.splatter_noise(full, rng, lib_scale=LIB_SCALE, bcv=BCV,
                                   dropout_mid=DROPOUT_MID)
            for p, _ in PIPELINES:
                x, y = preprocess(obs, p)
                meas = A.all_measures(x, y, n_kernel=1000, seed=r)
                for m in METHODS:
                    acc[(p, m)].append(meas[m])
        for p, _ in PIPELINES:
            for m in METHODS:
                v = np.array(acc[(p, m)], float)
                out[f'{label}__{p}__{m}__mean'] = np.array([np.nanmean(v)])
                out[f'{label}__{p}__{m}__std'] = np.array([np.nanstd(v)])
            print(f"[{label}] {p:9s} | " +
                  " ".join(f"{m}={out[f'{label}__{p}__{m}__mean'][0]:+.3f}" for m in METHODS))

        # ---- permutation z-scores per pipeline ----
        zacc = {(p, m): [] for p, _ in PIPELINES for m in METHODS}
        for r in range(N_REP_Z):
            rng = np.random.default_rng(9100 + 31 * r)
            full = np.concatenate([C, make_background(N_CELLS, rng)], axis=1)
            obs = A.splatter_noise(full, rng, lib_scale=LIB_SCALE, bcv=BCV,
                                   dropout_mid=DROPOUT_MID)
            for p, _ in PIPELINES:
                x, y = preprocess(obs, p)
                for m in METHODS:
                    _, z, _ = A.permutation_z(x, y, method=m, B=PERM_B,
                                              n_sub=PERM_NSUB, seed=1)
                    zacc[(p, m)].append(z)
        for p, _ in PIPELINES:
            for m in METHODS:
                out[f'{label}__{p}__{m}__z'] = np.array([np.mean(zacc[(p, m)])])
            print(f"[{label}] {p:9s} perm-z | " +
                  " ".join(f"{m}={out[f'{label}__{p}__{m}__z'][0]:+.1f}" for m in METHODS))
    return out


def plot(D):
    plotstyle.apply()
    keys = [p for p, _ in PIPELINES]
    labs = [l for _, l in PIPELINES]
    x = np.arange(len(keys)); w = 0.38

    # Layout: A-D as a single horizontal strip of four measure panels (like the
    # earlier Splatter figure), E-F as two wider z-score panels below.
    fig = plt.figure(figsize=(23, 10))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 1.05], hspace=0.42, wspace=0.32)
    panel = ['A', 'B', 'C', 'D']

    for k, m in enumerate(METHODS):
        ax = fig.add_subplot(gs[0, k])
        for j, (nut, label, col) in enumerate(CONDS):
            mu = [D[f'{label}__{p}__{m}__mean'][0] for p in keys]
            sd = [D[f'{label}__{p}__{m}__std'][0] for p in keys]
            ax.bar(x + (j - 0.5) * w, mu, w, yerr=sd, capsize=3, color=col,
                   label=label if k == 0 else None)
            ax.axhline(D[f'{label}__clean__{m}'][0], color=col, ls=':', lw=2.2,
                       label=(f'{label}: clean-count truth' if k == 0 else None))
        ax.axhline(0, color='k', lw=0.6)
        ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=9, rotation=30, ha='right')
        ax.set_ylabel(MLAB[m], fontsize=12)
        ax.set_title(f"({panel[k]}) {MLAB[m]}", fontsize=13)
        ax.grid(alpha=.3, axis='y')
        if k == 0:
            ax.legend(fontsize=9, loc='upper right')

    # (E) detection z in the truly-dependent condition (spans two columns)
    ax = fig.add_subplot(gs[1, 0:2])
    for m in METHODS:
        ax.plot(x, [D[f'CLOSED loop__{p}__{m}__z'][0] for p in keys], 'o-', lw=2, ms=8, label=MLAB[m])
    ax.axhline(2, color='green', ls='--', lw=1.2, label='z = 2  (~p < 0.05)')
    ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=10)
    ax.set_ylabel("permutation z-score")
    ax.set_title("(E) CLOSED loop - TRUE signal (higher is better)")
    ax.legend(fontsize=9, ncol=2); ax.grid(alpha=.3)

    # (F) false-positive z in the open loop (spans two columns)
    ax = fig.add_subplot(gs[1, 2:4])
    for m in METHODS:
        ax.plot(x, [D[f'NUTLIN-3__{p}__{m}__z'][0] for p in keys], 's-', lw=2, ms=8, label=MLAB[m])
    ax.axhline(2, color='green', ls='--', lw=1.2, label='z = 2  (~p < 0.05)')
    ax.axhline(0, color='k', lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=10)
    ax.set_ylabel("permutation z-score")
    ax.set_title("(F) NUTLIN-3 - library-size FALSE POSITIVE (lower is better)")
    ax.legend(fontsize=9, ncol=2); ax.grid(alpha=.3)

    fig.suptitle("Does standard preprocessing rescue the dependence measures under Splatter noise?\n"
                 f"[N={N_CELLS} cells + {N_BG} background genes, snapshot t={T_SNAP:.0f} min, "
                 f"dropout_mid={DROPOUT_MID:+.0f}; {N_REP} noise draws (A-D), {N_REP_Z} draws x "
                 f"{PERM_B} permutations, n_sub={PERM_NSUB} (E-F)]", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(str(FIG / "fig_step6b_normalization.png"), dpi=115)
    print("Saved fig_step6b_normalization.png")


def main():
    D = load_or_compute(DATA / "cache_step6b.npz", compute)
    plot(D)
    S = load_or_compute(DATA / "cache_step6b_ncells.npz", compute_ncells)
    plot_ncells(S)


if __name__ == "__main__":
    main()
