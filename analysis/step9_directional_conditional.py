# -*- coding: utf-8 -*-
"""
STEP 9 - Make Cor/MI/HSIC/dCor DIRECTIONAL (time) and separate regulation from correlation
(conditioning). Lagged Cor/dCor, Granger causality and partial correlation.
Results cached to data/cache_step9.npz; re-plotting only reloads.

Run:            python analysis/step9_directional_conditional.py
Force recompute: RECOMPUTE=1 python analysis/step9_directional_conditional.py
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
import directional as D
import plotstyle
from cache import load_or_compute

N_CELLS = 4000
HORIZON = 1000.0
DT = 2.0
LATE = 200.0
LAGS = [-12, -8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8, 12]
THREE_GENE_CACHE = DATA / "cache_step9_three_gene_nutlin.npz"


def direction(nutlin, seed):
    T = np.arange(0.0, HORIZON, DT)
    d = F.simulate_ensemble(N_CELLS, T, nutlin=nutlin, seed=seed)
    w = T >= LATE
    A = D.detrend_by_time(d[:, M.IDX['p53_mRNA'], :][:, w])
    B = D.detrend_by_time(d[:, M.IDX['Mdm2_mRNA'], :][:, w])
    return (np.array([D.lagged_xcorr(A, B, L) for L in LAGS]),
            np.array([D.lagged_xcorr(B, A, L) for L in LAGS]),
            np.array([D.lagged_dcor(A, B, L) for L in LAGS]),
            np.array([D.lagged_dcor(B, A, L) for L in LAGS]),
            np.array(D.granger_1step(A, B)), np.array(D.granger_1step(B, A)))


def common_driver(n=4000, seed=0):
    rng = np.random.default_rng(seed)
    Av = rng.gamma(3.0, 4.0, n); Bv = rng.poisson(0.8 * Av); Cv = rng.poisson(0.8 * Av)
    return np.corrcoef(Bv, Cv)[0, 1], D.partial_corr(Bv.astype(float), Cv.astype(float), Av)


def realdata():
    try:
        d = np.load(str(DATA / "data_step7_panel.npz"), allow_pickle=True)
    except FileNotFoundError:
        return np.nan, np.nan
    c = "Idasanutlin"; genes = list(d[f"{c}__genes"]); gi = {g: i for i, g in enumerate(genes)}
    X = d[f"{c}__X"].astype(float); tot = d[f"{c}__tot"].astype(float); cl = d[f"{c}__cline"]
    norm = np.log1p(X / tot[:, None] * 1e4); mask = cl == "LNCAPCLONEFGC_PROSTATE"
    col = lambda g: norm[mask, gi[g]]
    Z = np.mean([col(g) for g in ["FDXR", "GDF15", "SESN1", "DDB2", "TP53I3", "BTG2"]], axis=0)
    return np.corrcoef(col("MDM2"), col("CDKN1A"))[0, 1], D.partial_corr(col("MDM2"), col("CDKN1A"), Z)


def compute():
    out = {'lags_min': np.array(LAGS) * DT}
    for label, nut, seed in [("closed", 0.0, 91), ("nutlin", 1.0, 92)]:
        xc_ab, xc_ba, dc_ab, dc_ba, gab, gba = direction(nut, seed)
        out.update({f'{label}__xc_ab': xc_ab, f'{label}__xc_ba': xc_ba,
                    f'{label}__dc_ab': dc_ab, f'{label}__dc_ba': dc_ba,
                    f'{label}__gab': gab, f'{label}__gba': gba})
        print(f"[{label}] Granger A->B={float(gab):.3f} B->A={float(gba):.3f}")
    ms, ps = common_driver(); rm, rp = realdata()
    out['cd'] = np.array([ms, ps]); out['real'] = np.array([rm, rp])
    print(f"[common-driver] corr={ms:+.3f} partial={ps:+.3f} | [real] corr={rm:+.3f} partial={rp:+.3f}")
    return out


def plot(D_):
    plotstyle.apply()
    lags = D_['lags_min']
    fig = plt.figure(figsize=(13, 15))
    gs = fig.add_gridspec(3, 2)
    axs = np.array([[fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])],
                    [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]],
                   dtype=object)
    # Only tau >= 0 is needed: the two directional curves are DISTINCT there
    # (C_BA(tau) = C_AB(-tau)). At positive lag, A->B high while B->A low => A leads B.
    pos = lags >= 0
    lp = lags[pos]

    def _direction_panel(ax, ab, ba, ylabel, title):
        ab = np.asarray(ab)[pos]; ba = np.asarray(ba)[pos]
        ax.plot(lp, ab, 'o-', lw=2.6, color='C2', label='p53(t) -> MDM2(t+tau)')
        ax.plot(lp, ba, 's-', lw=2.2, color='C4', label='MDM2(t) -> p53(t+tau)')
        kpk = int(np.argmax(ab))
        ax.plot(lp[kpk], ab[kpk], '*', ms=18, color='crimson', zorder=5,
                label=f'p53->MDM2 peak at tau = {lp[kpk]:.0f} min')
        ax.set_xlabel('lag tau >= 0 (min)'); ax.set_ylabel(ylabel)
        ax.set_title(title); ax.legend(fontsize=10, loc='upper right'); ax.grid(alpha=.3)

    # Panel A: cross-correlation, with the negative-feedback overshoot annotated and the
    # Nutlin (open-loop) curve overlaid to show the lobe structure exists only when the loop is closed.
    axA = axs[0, 0]
    ab = np.asarray(D_['closed__xc_ab'])[pos]; ba = np.asarray(D_['closed__xc_ba'])[pos]
    abn = np.asarray(D_['nutlin__xc_ab'])[pos]
    axA.axhline(0, color='gray', lw=0.8)
    axA.plot(lp, ab, 'o-', lw=2.6, color='C2', label='p53(t)->MDM2(t+tau)  [closed]')
    axA.plot(lp, ba, 's-', lw=2.0, color='C4', label='MDM2(t)->p53(t+tau)  [closed]')
    axA.plot(lp, abn, '--', lw=1.8, color='gray', label='p53->MDM2  [Nutlin, open loop]')
    kpk = int(np.argmax(ab)); kmin = int(np.argmin(ab))
    axA.plot(lp[kpk], ab[kpk], '*', ms=18, color='crimson', zorder=5, label=f'peak at tau = {lp[kpk]:.0f} min')
    # Keep the label directly above the minimum: the arrow is vertical and lands
    # exactly on the negative-lobe point, avoiding the earlier diagonal callout.
    axA.annotate('negative lobe =\nnegative-feedback overshoot\n(closed loop only)',
                 xy=(lp[kmin], ab[kmin]), xytext=(lp[kmin], 0.22),
                 fontsize=9, color='C3', ha='center', va='bottom',
                 arrowprops=dict(arrowstyle='->', color='C3', lw=1.5,
                                 connectionstyle='arc3,rad=0'))
    axA.set_xlabel('lag tau >= 0 (min)'); axA.set_ylabel('cross-correlation (fluctuations)')
    axA.set_title("(A) Lagged cross-correlation (closed loop)\np53 leads MDM2; negative lobe = feedback overshoot")
    axA.legend(fontsize=9, loc='upper right'); axA.grid(alpha=.3)

    _direction_panel(axs[0, 1], D_['closed__dc_ab'], D_['closed__dc_ba'], 'lagged distance correlation',
                     "(B) Lagged distance correlation (closed loop, nonlinear)\ndCor>=0 by design: cannot show the negative lobe")
    ax = axs[1, 0]
    x = np.arange(2); w = 0.32
    forward = [float(D_['closed__gab']), float(D_['nutlin__gab'])]
    reverse = [float(D_['closed__gba']), float(D_['nutlin__gba'])]
    bars_f = ax.bar(x-w/2, forward, w, color='C2', label='Granger p53->MDM2')
    bars_r = ax.bar(x+w/2, reverse, w, color='C4', label='Granger MDM2->p53')
    ax.set_xticks(x); ax.set_xticklabels(['CLOSED', 'NUTLIN'])
    ax.set_ylabel('Granger strength  log(var reduced / var full)')
    ax.set_ylim(0, max(forward + reverse) * 1.25)
    ax.legend(fontsize=10)
    ax.set_title("(C) One-step Granger causality:\np53->MDM2 is recovered in the dynamic closed loop")
    def granger_label(value):
        # Values below 0.001 are real but would be hidden by fixed 3-decimal rounding.
        return f'{value:.3f}' if abs(value) >= 1e-3 else f'{value:.2e}'

    for bars in (bars_f, bars_r):
        for bar in bars:
            value = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, value + 0.006,
                    granger_label(value), ha='center', va='bottom', fontsize=9)
    ax.grid(alpha=.3, axis='y')
    rm, rp = D_['real']
    ax = axs[1, 1]
    t3 = D_['fork3_times']; logfc3 = D_['fork3_mean_logfc']
    ax.axhline(0, color='gray', lw=0.8)
    ax.plot(t3, logfc3[0], 'o-', lw=2.6, color='C3', label='p53 protein')
    ax.plot(t3, logfc3[1], 'o-', lw=2.3, color='C4', label='TP53 mRNA')
    ax.plot(t3, logfc3[2], 'o-', lw=2.6, color='C0', label='MDM2 mRNA')
    ax.plot(t3, logfc3[3], 'o-', lw=2.6, color='C2', label='CDKN1A mRNA')
    ax.set_xlabel('time after Nutlin (min)')
    ax.set_ylabel('log2 fold-change from t = 0')
    ax.set_title("(D) Open-loop response under Nutlin:\n"
                 "TP53 mRNA stays flat; p53 protein and its targets rise")
    ax.legend(fontsize=8)
    ax.grid(alpha=.3, axis='y')
    # Panel E: synthetic clean/noisy and real-data reference. Each synthetic cell is
    # observed once at a different response phase; p53 is the true protein-level
    # common driver in this six-species simulation. The Splatter version adds
    # library-size variation, overdispersion and expression-dependent dropout to
    # the observed target transcripts.
    ax = fig.add_subplot(gs[2, :])
    synthetic = np.r_[D_['fork3_snapshot_stats'],
                      D_['fork3_snapshot_splatter_stats'],
                      D_['real']]
    xpos = np.arange(6)
    ax.bar(xpos, synthetic, color=['C0', 'C3', 'C0', 'C3', 'C0', 'C3'])
    ax.axhline(0, color='gray', lw=0.8)
    ax.set_xticks(xpos)
    ax.set_xticklabels(['clean\ncorr',
                        'clean\npartial | p53',
                        'Splatter\ncorr',
                        'Splatter\npartial | p53',
                        'real\ncorr',
                        'real\npartial | p53-proxy'])
    ax.set_ylabel('correlation')
    ax.set_ylim(-0.08, 0.8)
    if 'fork3_snapshot_splatter_zero_rates' in D_:
        zero_note = (f"Splatter zeros: MDM2 {100 * D_['fork3_snapshot_splatter_zero_rates'][0]:.1f}%, "
                     f"CDKN1A {100 * D_['fork3_snapshot_splatter_zero_rates'][1]:.1f}%")
    else:
        zero_note = "technical noise leaves residual association"
    ax.set_title("(E) Marginal vs partial correlation:\n"
                 f"clean synthetic, Splatter-noised synthetic, and real idasanutlin data ({zero_note})")
    ax.grid(alpha=.3, axis='y')
    for xj, value in zip(xpos, synthetic):
        va = 'bottom' if value >= 0 else 'top'
        offset = 0.02 if value >= 0 else -0.02
        ax.text(xj, value + offset, f'{value:+.3f}', ha='center', va=va, fontsize=10)
    fig.suptitle("Step 9 - Directional and conditional dependence for the p53-MDM2 mRNA pair\n"
                 f"[N={N_CELLS} cells per condition]", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(str(FIG / "fig_step9_directional.png"), dpi=120)
    print("Saved fig_step9_directional.png")


def main():
    cache_path = DATA / "cache_step9.npz"
    D_ = load_or_compute(cache_path, compute)
    if not THREE_GENE_CACHE.exists():
        print("Missing three-gene cache - run analysis/step9b_three_gene_nutlin.py first.")
        return
    with np.load(THREE_GENE_CACHE, allow_pickle=True) as s3:
        D_['fork3_times'] = s3['times']
        species = {str(name): i for i, name in enumerate(s3['species'])}
        trajectories = s3['trajectories']
        means = np.array([
            trajectories[:, species['p53'], :].mean(axis=0),
            trajectories[:, species['p53_mRNA'], :].mean(axis=0),
            trajectories[:, species['Mdm2_mRNA'], :].mean(axis=0),
            trajectories[:, species['CDKN1A_mRNA'], :].mean(axis=0),
        ])
        D_['fork3_mean_logfc'] = np.log2((means + 1.0) / (means[:, :1] + 1.0))
        D_['fork3_snapshot_stats'] = s3['snapshot_stats']
        if 'snapshot_splatter_stats' not in s3:
            print("Missing Splatter snapshot fields - run analysis/step9b_three_gene_nutlin.py first.")
            return
        D_['fork3_snapshot_splatter_stats'] = s3['snapshot_splatter_stats']
        if 'snapshot_splatter_zero_rates' in s3:
            D_['fork3_snapshot_splatter_zero_rates'] = s3['snapshot_splatter_zero_rates']
    plot(D_)


if __name__ == "__main__":
    main()
