# -*- coding: utf-8 -*-
"""
STEP 7d - The DIRECT (TP53, MDM2) pair on real scRNA-seq, shown honestly.

Our ideal observable is the loop pair (p53, MDM2). On real scRNA-seq, however, TP53 mRNA is a
poor readout of p53 activity (p53 is regulated post-translationally), so the (TP53 mRNA,
MDM2 mRNA) dependence is weak. This step shows that directly and contrasts it with the p53-target
pair (MDM2, CDKN1A), motivating the choice made in step 7b.

Uses data/data_step7_panel.npz (6 h). WT line = LNCaP.
Run:  python analysis/step7d_direct_pair_realdata.py
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

WT_LINE = "LNCAPCLONEFGC_PROSTATE"
CONDS = [("DMSO", 'C0'), ("Idasanutlin", 'C3')]


def load():
    d = np.load(str(DATA / "data_step7_panel.npz"), allow_pickle=True)
    out = {}
    for c, _ in CONDS:
        genes = list(d[f"{c}__genes"]); X = d[f"{c}__X"].astype(float)
        tot = d[f"{c}__tot"].astype(float); cline = d[f"{c}__cline"]
        norm = np.log1p(X / tot[:, None] * 1e4)
        out[c] = dict(gidx={g: i for i, g in enumerate(genes)}, norm=norm, cline=cline)
    return out


def expr(panel, cond, gene):
    dd = panel[cond]; mask = dd["cline"] == WT_LINE
    return dd["norm"][mask, dd["gidx"][gene]]


def main():
    p = load()
    plotstyle.apply()
    fig, axs = plt.subplots(1, 3, figsize=(19, 5.5))

    # (A) mean expression: TP53 mRNA barely changes; MDM2/CDKN1A rise
    genes = ["TP53", "MDM2", "CDKN1A"]
    x = np.arange(len(genes)); w = 0.38
    axs[0].bar(x - w/2, [expr(p, "DMSO", g).mean() for g in genes], w, label="DMSO", color='C0')
    axs[0].bar(x + w/2, [expr(p, "Idasanutlin", g).mean() for g in genes], w, label="Idasanutlin", color='C3')
    axs[0].set_xticks(x); axs[0].set_xticklabels(genes)
    axs[0].set_ylabel("mean expression (log-norm)"); axs[0].legend()
    axs[0].set_title("(A) TP53 mRNA hardly moves (post-translational);\nMDM2 & CDKN1A are induced")

    # (B) scatter TP53 vs MDM2 (the direct loop pair) - weak/no structure
    for cond, col in CONDS:
        xx = expr(p, cond, "TP53"); yy = expr(p, cond, "MDM2")
        axs[1].scatter(xx, yy, s=10, alpha=0.4, color=col,
                       label=f"{cond[:5]} (dCor={A.distance_corr(xx, yy):.2f})")
    axs[1].set_xlabel("TP53 mRNA (log-norm)"); axs[1].set_ylabel("MDM2 mRNA (log-norm)")
    axs[1].set_title("(B) Direct loop pair TP53-MDM2:\nweak dependence (TP53 mRNA uninformative)"); axs[1].legend()

    # (C) dCor: TP53-MDM2 vs MDM2-CDKN1A, DMSO vs Idasanutlin
    pairs = [("TP53", "MDM2"), ("MDM2", "CDKN1A")]
    xg = np.arange(len(pairs))
    dmso = [A.distance_corr(expr(p, "DMSO", a), expr(p, "DMSO", b)) for a, b in pairs]
    idasa = [A.distance_corr(expr(p, "Idasanutlin", a), expr(p, "Idasanutlin", b)) for a, b in pairs]
    axs[2].bar(xg - w/2, dmso, w, label="DMSO", color='C0')
    axs[2].bar(xg + w/2, idasa, w, label="Idasanutlin", color='C3')
    axs[2].set_xticks(xg); axs[2].set_xticklabels([f"{a}-{b}" for a, b in pairs])
    axs[2].set_ylabel("dCor"); axs[2].legend()
    axs[2].set_title("(C) Direct pair stays weak; the p53-target\npair MDM2-CDKN1A carries the signal")

    n = int(np.sum(p["Idasanutlin"]["cline"] == WT_LINE))
    fig.suptitle(f"Step 7d - Why not the direct (TP53, MDM2) pair on real data? WT line {WT_LINE.split('_')[0]} "
                 f"(6 h, n={n} treated)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(str(FIG / "fig_step7d_direct_pair.png"), dpi=120)
    print("TP53-MDM2 dCor:", {c: round(A.distance_corr(expr(p, c, 'TP53'), expr(p, c, 'MDM2')), 3) for c, _ in CONDS})
    print("MDM2-CDKN1A dCor:", {c: round(A.distance_corr(expr(p, c, 'MDM2'), expr(p, c, 'CDKN1A')), 3) for c, _ in CONDS})
    print("Saved fig_step7d_direct_pair.png")


if __name__ == "__main__":
    main()
