# -*- coding: utf-8 -*-
"""
STEP 7c - 6 h vs 24 h timepoint comparison for idasanutlin (MIX-seq, McFarland et al. 2020).

Idasanutlin (2.5 uM) was assayed at 6 h and 24 h. This step asks which timepoint shows the
stronger p53 response, in the TP53 wild-type line, by comparing:
  (A) induction of p53 target genes (Idasanutlin - DMSO mean, log-norm) at 6 h vs 24 h;
  (B) the MDM2-CDKN1A co-expression (dCor) at 6 h vs 24 h.

Needs data/data_step7_panel.npz (6 h) and data/data_step7_panel_24hr.npz (24 h) from step 7a.
Run:  python analysis/step7c_timepoint_comparison.py
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

WT_LINE = "LNCAPCLONEFGC_PROSTATE"   # top TP53-WT responder (see step 7b ranking)
TARGETS = ["MDM2", "CDKN1A", "FDXR", "TP53I3", "BTG2", "GDF15", "SESN1", "DDB2"]
TPS = [("6hr", "data_step7_panel.npz"), ("24hr", "data_step7_panel_24hr.npz")]


def load(fname):
    d = np.load(str(DATA / fname), allow_pickle=True)
    out = {}
    for c in ("DMSO", "Idasanutlin"):
        genes = list(d[f"{c}__genes"]); X = d[f"{c}__X"].astype(float)
        tot = d[f"{c}__tot"].astype(float); cline = d[f"{c}__cline"]
        norm = np.log1p(X / tot[:, None] * 1e4)
        out[c] = dict(gidx={g: i for i, g in enumerate(genes)}, norm=norm, cline=cline)
    return out


def expr(panel, cond, gene):
    dd = panel[cond]; mask = dd["cline"] == WT_LINE
    return dd["norm"][mask, dd["gidx"][gene]]


def main():
    panels = {}
    for tp, fname in TPS:
        try:
            panels[tp] = load(fname)
        except FileNotFoundError:
            print(f"Missing {fname} - run step7a first."); return

    # (A) induction per target gene at each timepoint
    induction = {}
    ncells = {}
    for tp, _ in TPS:
        p = panels[tp]
        induction[tp] = [expr(p, "Idasanutlin", g).mean() - expr(p, "DMSO", g).mean() for g in TARGETS]
        ncells[tp] = (int(np.sum(p["DMSO"]["cline"] == WT_LINE)),
                      int(np.sum(p["Idasanutlin"]["cline"] == WT_LINE)))
    # (B) MDM2-CDKN1A dCor
    dcor = {}
    for tp, _ in TPS:
        p = panels[tp]
        dcor[tp] = {}
        for cond in ("DMSO", "Idasanutlin"):
            dcor[tp][cond] = A.distance_corr(expr(p, cond, "MDM2"), expr(p, cond, "CDKN1A"))

    print(f"WT line {WT_LINE}; cells (DMSO/Idasa): 6h={ncells['6hr']}, 24h={ncells['24hr']}")
    print("Mean induction (Idasa-DMSO) summed over targets:",
          {tp: round(float(np.sum(induction[tp])), 2) for tp, _ in TPS})
    print("MDM2-CDKN1A dCor:", {tp: {k: round(v, 2) for k, v in dcor[tp].items()} for tp, _ in TPS})

    plotstyle.apply()
    fig, axs = plt.subplots(1, 2, figsize=(15, 6))
    x = np.arange(len(TARGETS)); w = 0.38
    axs[0].bar(x - w/2, induction["6hr"], w, label="6 h", color='C0')
    axs[0].bar(x + w/2, induction["24hr"], w, label="24 h", color='C1')
    axs[0].set_xticks(x); axs[0].set_xticklabels(TARGETS, rotation=45, ha='right')
    axs[0].set_ylabel("induction  (Idasanutlin - DMSO, log-norm)")
    axs[0].set_title("(A) p53-target induction: stronger at 24 h"); axs[0].legend(); axs[0].grid(alpha=.3, axis='y')

    x2 = np.arange(2);
    axs[1].bar(x2 - w/2, [dcor["6hr"]["DMSO"], dcor["24hr"]["DMSO"]], w, label="DMSO", color='C0')
    axs[1].bar(x2 + w/2, [dcor["6hr"]["Idasanutlin"], dcor["24hr"]["Idasanutlin"]], w, label="Idasanutlin", color='C3')
    axs[1].set_xticks(x2); axs[1].set_xticklabels(["6 h", "24 h"])
    axs[1].set_ylabel("dCor (MDM2, CDKN1A)")
    axs[1].set_title("(B) MDM2-CDKN1A co-expression: 6 h vs 24 h"); axs[1].legend(); axs[1].grid(alpha=.3, axis='y')

    fig.suptitle(f"Step 7c - Idasanutlin (2.5 uM) timepoint comparison, WT line {WT_LINE.split('_')[0]} "
                 f"(6 h: n={ncells['6hr'][1]} treated; 24 h: n={ncells['24hr'][1]} treated)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(str(FIG / "fig_step7c_timepoints.png"), dpi=120)
    print("Saved fig_step7c_timepoints.png")


if __name__ == "__main__":
    main()
