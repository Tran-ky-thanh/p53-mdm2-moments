# -*- coding: utf-8 -*-
"""
STEP 7e - Generalize the 2x2 control (Section 4.7 / Figure 8E) across the whole cell-line panel.

For every cell line with enough QC-passing cells in BOTH conditions (6 h), measure:
  - the p53 response = induction of MDM2 + CDKN1A (Idasanutlin - DMSO, log-norm mean), and
  - the MDM2-CDKN1A co-expression (distance correlation) in DMSO vs Idasanutlin.
Question: do the responsive (TP53-WT-like) lines gain MDM2-CDKN1A co-expression under the drug?

Uses data/data_step7_panel.npz (6 h) from step 7a.
Run:  python analysis/step7e_across_cell_lines.py
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
FIG = ROOT / "figures"; DATA = ROOT / "data"

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import association as A
import plotstyle

MIN_CELLS = 40
RESP_HI = 1.0    # p53 response threshold for "responder"
RESP_LO = 0.5    # below this = "non-responder"


def load(cond):
    d = np.load(str(DATA / "data_step7_panel.npz"), allow_pickle=True)
    g = list(d[f"{cond}__genes"]); gi = {x: i for i, x in enumerate(g)}
    X = d[f"{cond}__X"].astype(float); tot = d[f"{cond}__tot"].astype(float); cl = d[f"{cond}__cline"]
    return gi, np.log1p(X / tot[:, None] * 1e4), cl


def main():
    giD, nD, clD = load("DMSO"); giI, nI, clI = load("Idasanutlin")
    col = lambda gi, n, cl, L, g: n[cl == L, gi[g]]
    rows = []
    for L in set(clD) & set(clI):
        nd = int((clD == L).sum()); ni = int((clI == L).sum())
        if nd < MIN_CELLS or ni < MIN_CELLS:
            continue
        resp = sum(col(giI, nI, clI, L, g).mean() - col(giD, nD, clD, L, g).mean()
                   for g in ["MDM2", "CDKN1A"])
        dc_d = A.distance_corr(col(giD, nD, clD, L, "MDM2"), col(giD, nD, clD, L, "CDKN1A"))
        dc_i = A.distance_corr(col(giI, nI, clI, L, "MDM2"), col(giI, nI, clI, L, "CDKN1A"))
        rows.append(dict(line=L.split('_')[0], resp=resp, dc_d=dc_d, dc_i=dc_i, gain=dc_i - dc_d))
    rows.sort(key=lambda r: -r['resp'])
    resp = np.array([r['resp'] for r in rows]); gain = np.array([r['gain'] for r in rows])
    rr = np.corrcoef(resp, gain)[0, 1]
    print(f"{len(rows)} cell lines; corr(p53 response, dCor gain) = {rr:+.2f}")
    for r in rows:
        print(f"  {r['line']:24s} resp={r['resp']:+.2f}  dCor {r['dc_d']:+.2f}->{r['dc_i']:+.2f}  gain={r['gain']:+.2f}")

    plotstyle.apply()
    fig, ax = plt.subplots(figsize=(9, 7.5))

    # Standalone panel: p53 response vs dCor gain, one point per line
    ax.scatter(resp, gain, s=70, color='C0', zorder=3)
    for r in rows:
        if r['resp'] > RESP_HI or abs(r['gain']) > 0.2:
            ax.annotate(r['line'], (r['resp'], r['gain']), fontsize=8,
                        xytext=(4, 4), textcoords='offset points')
    b, a = np.polyfit(resp, gain, 1)
    xs = np.linspace(resp.min(), resp.max(), 50)
    ax.plot(xs, a + b * xs, '--', color='C3', label=f'trend (r = {rr:+.2f})')
    ax.axhline(0, color='gray', lw=0.6)
    ax.set_xlabel("p53 response  (MDM2+CDKN1A induction, Idasanutlin - DMSO)")
    ax.set_ylabel("dCor gain  (Idasanutlin - DMSO)")
    ax.set_title(f"Step 7e - Across {len(rows)} MIX-seq cell lines (6 h, >= {MIN_CELLS} cells/condition):\n"
                 "more p53 response -> more MDM2-CDKN1A co-expression gain")
    ax.legend(); ax.grid(alpha=.3)

    fig.tight_layout()
    fig.savefig(str(FIG / "fig_step7e_across_lines.png"), dpi=120)
    print("Saved fig_step7e_across_lines.png")


if __name__ == "__main__":
    main()
