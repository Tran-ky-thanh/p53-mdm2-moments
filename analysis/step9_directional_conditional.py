# -*- coding: utf-8 -*-
"""
STEP 9 - Make Cor/MI/HSIC/dCor DIRECTIONAL (time) and separate regulation from correlation
(conditioning). Lagged Cor/dCor, Granger, Transfer Entropy; partial correlation.
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
            np.array(D.granger_1step(A, B)), np.array(D.granger_1step(B, A)),
            np.array(D.transfer_entropy(A, B)), np.array(D.transfer_entropy(B, A)))


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
        xc_ab, xc_ba, dc_ab, dc_ba, gab, gba, teab, teba = direction(nut, seed)
        out.update({f'{label}__xc_ab': xc_ab, f'{label}__xc_ba': xc_ba,
                    f'{label}__dc_ab': dc_ab, f'{label}__dc_ba': dc_ba,
                    f'{label}__gab': gab, f'{label}__gba': gba, f'{label}__teab': teab, f'{label}__teba': teba})
        print(f"[{label}] Granger A->B={float(gab):.3f} B->A={float(gba):.3f} | TE A->B={float(teab):.4f} B->A={float(teba):.4f}")
    ms, ps = common_driver(); rm, rp = realdata()
    out['cd'] = np.array([ms, ps]); out['real'] = np.array([rm, rp])
    print(f"[common-driver] corr={ms:+.3f} partial={ps:+.3f} | [real] corr={rm:+.3f} partial={rp:+.3f}")
    return out


def plot(D_):
    plotstyle.apply()
    lags = D_['lags_min']
    fig, axs = plt.subplots(3, 2, figsize=(13, 15))
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
    axA.annotate('negative lobe =\nnegative-feedback overshoot\n(closed loop only)',
                 xy=(lp[kmin], ab[kmin]), xytext=(lp[kmin]-9, 0.22), fontsize=9, color='C3', ha='left',
                 arrowprops=dict(arrowstyle='->', color='C3', lw=1.5))
    axA.set_xlabel('lag tau >= 0 (min)'); axA.set_ylabel('cross-correlation (fluctuations)')
    axA.set_title("(A) Lagged cross-correlation (closed loop)\np53 leads MDM2; negative lobe = feedback overshoot")
    axA.legend(fontsize=9, loc='upper right'); axA.grid(alpha=.3)

    _direction_panel(axs[0, 1], D_['closed__dc_ab'], D_['closed__dc_ba'], 'lagged distance correlation',
                     "(B) Lagged distance correlation (closed loop, nonlinear)\ndCor>=0 by design: cannot show the negative lobe")
    ax = axs[1, 0]
    x = np.arange(2); w = 0.2
    ax.bar(x-1.5*w, [float(D_['closed__gab']), float(D_['nutlin__gab'])], w, color='C2', label='Granger A->B')
    ax.bar(x-0.5*w, [float(D_['closed__gba']), float(D_['nutlin__gba'])], w, color='C4', label='Granger B->A')
    ax.bar(x+0.5*w, [float(D_['closed__teab'])*10, float(D_['nutlin__teab'])*10], w, color='C2', alpha=0.5, label='TE A->B x10')
    ax.bar(x+1.5*w, [float(D_['closed__teba'])*10, float(D_['nutlin__teba'])*10], w, color='C4', alpha=0.5, label='TE B->A x10')
    ax.set_xticks(x); ax.set_xticklabels(['CLOSED', 'NUTLIN'])
    ax.set_ylabel('directional strength'); ax.legend(fontsize=10)
    ax.set_title("(C) Granger & transfer entropy:\nA->B >> B->A  (direction p53->MDM2)")
    ax.grid(alpha=.3, axis='y')
    ms, ps = D_['cd']; rm, rp = D_['real']
    ax = axs[1, 1]
    ax.bar([0, 1], [ms, ps], color=['C0', 'C3'])
    ax.set_xticks([0, 1]); ax.set_xticklabels(['corr(B,C)\nmarginal', 'partial(B,C | A)'])
    ax.set_ylabel('correlation')
    ax.set_title("(D) Common driver B<-A->C:\npartial ~ 0 => not a direct edge")
    ax.grid(alpha=.3, axis='y')
    ax = axs[2, 0]
    ax.bar([0, 1], [rm, rp], color=['C0', 'C3'])
    ax.set_xticks([0, 1]); ax.set_xticklabels(['corr(MDM2,CDKN1A)', 'partial(.|p53-act)'])
    ax.set_ylabel('correlation')
    ax.set_title("(E) Real data: partial correlation barely drops\n(imperfect mRNA proxy of p53 => inconclusive)")
    ax.grid(alpha=.3, axis='y')
    axs[2, 1].axis('off')
    fig.suptitle("Step 9 - Directional and conditional dependence for the p53-MDM2 mRNA pair\n"
                 f"[N={N_CELLS} cells per condition]", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(str(FIG / "fig_step9_directional.png"), dpi=120)
    print("Saved fig_step9_directional.png")


def main():
    D_ = load_or_compute(DATA / "cache_step9.npz", compute)
    plot(D_)


if __name__ == "__main__":
    main()
