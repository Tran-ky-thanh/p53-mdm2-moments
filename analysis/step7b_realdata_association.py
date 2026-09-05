# -*- coding: utf-8 -*-
"""
STEP 7b - Apply our dependence measures to REAL data (MIX-seq) to test whether they
actually detect p53-MDM2 regulation, with negative controls from other gene pairs.

Uses data/data_step7_panel.npz from step 7a. Pipeline:
  1) Normalize CP10k + log1p (removes the library-size confounder highlighted in step 6).
  2) Identify TP53-WT lines: those where MDM2 & CDKN1A rise strongly under Idasanutlin
     (Nutlin activates p53 only in TP53-WT). Pick one WT line + one mutant line (control).
  3) Validate p53 activation: p53-target means rise in the WT line under Idasanutlin.
  4) Measure dependence (Pearson/Spearman/MI/HSIC/dCor) for:
       - p53-target pairs (co-regulated)  - expected dependence
       - the loop pair TP53-MDM2
       - NEGATIVE-CONTROL pairs (housekeeping/random) - expected ~0
     Compare DMSO (closed) vs Idasanutlin (open); 2x2 WT/mutant control; permutation z-scores.

Run:  python analysis/step7b_realdata_association.py
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
CONDS = ["DMSO", "Idasanutlin"]
TARGET_PAIRS = [("MDM2", "CDKN1A"), ("MDM2", "FDXR"), ("CDKN1A", "FDXR"),
                ("FDXR", "BTG2"), ("MDM2", "TP53")]
CONTROL_PAIRS = [("ACTB", "GAPDH"), ("B2M", "RPL13A"), ("GAPDH", "PGK1"),
                 ("ACTB", "B2M"), ("PPIA", "EEF1A1")]
METHODS = ['pearson', 'spearman', 'MI', 'HSIC', 'dCor']


def load():
    d = np.load(str(DATA / "data_step7_panel.npz"), allow_pickle=True)
    data = {}
    for c in CONDS:
        genes = list(d[f"{c}__genes"]); X = d[f"{c}__X"].astype(float)
        tot = d[f"{c}__tot"].astype(float); cline = d[f"{c}__cline"]
        norm = np.log1p(X / tot[:, None] * 1e4)
        data[c] = dict(genes=genes, gidx={g: i for i, g in enumerate(genes)}, norm=norm, cline=cline)
    return data


def find_lines(data):
    common = set(data["DMSO"]["cline"]) & set(data["Idasanutlin"]["cline"])
    rows = []
    for line in common:
        n_d = int(np.sum(data["DMSO"]["cline"] == line)); n_i = int(np.sum(data["Idasanutlin"]["cline"] == line))
        if n_d < MIN_CELLS or n_i < MIN_CELLS:
            continue
        delta = 0.0
        for g in ["MDM2", "CDKN1A"]:
            md = data["DMSO"]["norm"][data["DMSO"]["cline"] == line, data["DMSO"]["gidx"][g]].mean()
            mi = data["Idasanutlin"]["norm"][data["Idasanutlin"]["cline"] == line, data["Idasanutlin"]["gidx"][g]].mean()
            delta += mi - md
        rows.append((line, n_d, n_i, delta))
    rows.sort(key=lambda r: -r[3])
    return rows


def get_expr(data, cond, line, gene):
    dd = data[cond]; mask = dd["cline"] == line
    return dd["norm"][mask, dd["gidx"][gene]]


def measure_pairs(data, line, pairs):
    res = {}
    for (g1, g2) in pairs:
        res[(g1, g2)] = {}
        for cond in CONDS:
            x = get_expr(data, cond, line, g1); y = get_expr(data, cond, line, g2)
            res[(g1, g2)][cond] = {m: A.permutation_z(x, y, method=m, B=200,
                                                      n_sub=min(1500, len(x)), seed=0) for m in METHODS}
    return res


def main():
    data = load()
    ranking = find_lines(data)
    print("Cell lines ranked by p53 response (delta mean MDM2+CDKN1A):")
    for line, nd, ni, delta in ranking[:15]:
        print(f"  {line:32s} n_DMSO={nd:4d} n_Idasa={ni:4d}  delta={delta:+.3f}")
    wt_line = ranking[0][0]; mut_line = ranking[-1][0]
    print(f"\nWT (functional loop) = {wt_line}; mutant (control) = {mut_line}")

    targets = [g for g in ["MDM2", "CDKN1A", "FDXR", "TP53I3", "BTG2", "GDF15", "SESN1", "DDB2"]
               if g in data["DMSO"]["gidx"]]
    means = {}
    for line, tag in [(wt_line, "WT"), (mut_line, "MUT")]:
        for cond in CONDS:
            means[(tag, cond)] = [get_expr(data, cond, line, g).mean() for g in targets]

    res_t = measure_pairs(data, wt_line, TARGET_PAIRS)
    res_c = measure_pairs(data, wt_line, CONTROL_PAIRS)
    res_star = {}
    for line, tag in [(wt_line, "WT"), (mut_line, "MUT")]:
        for cond in CONDS:
            x = get_expr(data, cond, line, "MDM2"); y = get_expr(data, cond, line, "CDKN1A")
            res_star[(tag, cond)] = {m: A.permutation_z(x, y, method=m, B=300,
                                                        n_sub=min(1500, len(x)), seed=0) for m in METHODS}
    print("\n[MDM2-CDKN1A] 2x2 control (dCor, z):")
    for tag in ["WT", "MUT"]:
        for cond in CONDS:
            r = res_star[(tag, cond)]['dCor']
            print(f"  {tag} * {cond:12s}: dCor={r[0]:+.2f} z={r[1]:+.1f} p={r[2]:.3f}")

    plotstyle.apply()
    fig = plt.figure(figsize=(14, 17))
    ax = fig.add_subplot(3, 2, 1)
    x = np.arange(len(targets)); w = 0.2
    ax.bar(x - 1.5*w, means[("WT", "DMSO")], w, label="WT/DMSO", color='C0')
    ax.bar(x - 0.5*w, means[("WT", "Idasanutlin")], w, label="WT/Idasa", color='C3')
    ax.bar(x + 0.5*w, means[("MUT", "DMSO")], w, label="MUT/DMSO", color='C0', alpha=0.4)
    ax.bar(x + 1.5*w, means[("MUT", "Idasanutlin")], w, label="MUT/Idasa", color='C3', alpha=0.4)
    ax.set_xticks(x); ax.set_xticklabels(targets, rotation=45, ha='right', fontsize=11)
    ax.set_ylabel("mean log-norm expr"); ax.legend(fontsize=10)
    ax.set_title(f"(A) p53 activation: targets rise in WT ({wt_line.split('_')[0]});\nmutant unchanged", fontsize=13)

    ax = fig.add_subplot(3, 2, 2)
    for cond, col in [("DMSO", 'C0'), ("Idasanutlin", 'C3')]:
        x1 = get_expr(data, cond, wt_line, "MDM2"); y1 = get_expr(data, cond, wt_line, "CDKN1A")
        ax.scatter(x1, y1, s=8, alpha=0.4, color=col, label=f"{cond[:5]} (dCor={A.distance_corr(x1,y1):.2f})")
    ax.set_xlabel("MDM2 (log-norm)"); ax.set_ylabel("CDKN1A (log-norm)")
    ax.set_title("(B) MDM2-CDKN1A in WT: co-expression\nappears when the loop opens", fontsize=13); ax.legend(fontsize=11)

    ax = fig.add_subplot(3, 2, 3)
    for cond, col in [("DMSO", 'C0'), ("Idasanutlin", 'C3')]:
        x1 = get_expr(data, cond, mut_line, "MDM2"); y1 = get_expr(data, cond, mut_line, "CDKN1A")
        ax.scatter(x1, y1, s=8, alpha=0.4, color=col, label=f"{cond[:5]} (dCor={A.distance_corr(x1,y1):.2f})")
    ax.set_xlabel("MDM2 (log-norm)"); ax.set_ylabel("CDKN1A (log-norm)")
    ax.set_title(f"(C) MDM2-CDKN1A in MUTANT ({mut_line.split('_')[0]}):\nno co-expression (control)", fontsize=13)
    ax.legend(fontsize=11)

    ax = fig.add_subplot(3, 2, 4)
    xm = np.arange(len(METHODS)); w = 0.38
    ax.bar(xm - w/2, [res_star[("WT", "DMSO")][m][0] for m in METHODS], w, label="DMSO", color='C0')
    ax.bar(xm + w/2, [res_star[("WT", "Idasanutlin")][m][0] for m in METHODS], w, label="Idasa", color='C3')
    ax.set_xticks(xm); ax.set_xticklabels(['Pear', 'Spear', 'MI', 'HSIC', 'dCor'], fontsize=12)
    ax.set_ylabel("statistic"); ax.legend(fontsize=11)
    ax.set_title("(D) MDM2-CDKN1A (WT): every method\ndetects co-regulation when loop opens", fontsize=13)

    ax = fig.add_subplot(3, 2, 5)
    x2 = np.arange(2)
    ax.bar(x2 - w/2, [res_star[("WT", "DMSO")]['dCor'][0], res_star[("MUT", "DMSO")]['dCor'][0]], w,
           label="DMSO", color='C0')
    ax.bar(x2 + w/2, [res_star[("WT", "Idasanutlin")]['dCor'][0], res_star[("MUT", "Idasanutlin")]['dCor'][0]], w,
           label="Idasa", color='C3')
    ax.set_xticks(x2); ax.set_xticklabels([f"WT\n{wt_line.split('_')[0]}", f"MUT\n{mut_line.split('_')[0]}"], fontsize=12)
    ax.set_ylabel("dCor (MDM2-CDKN1A)"); ax.legend(fontsize=11)
    ax.set_title("(E) 2x2 control: only WT+Idasa\nshows strong co-regulation", fontsize=13)

    ax = fig.add_subplot(3, 2, 6)
    allpairs = TARGET_PAIRS + CONTROL_PAIRS
    delta = ([res_t[p]["Idasanutlin"]['dCor'][0] - res_t[p]["DMSO"]['dCor'][0] for p in TARGET_PAIRS] +
             [res_c[p]["Idasanutlin"]['dCor'][0] - res_c[p]["DMSO"]['dCor'][0] for p in CONTROL_PAIRS])
    cols = ['C3'] * len(TARGET_PAIRS) + ['gray'] * len(CONTROL_PAIRS)
    xx = np.arange(len(allpairs))
    ax.bar(xx, delta, color=cols); ax.axhline(0, color='k', lw=0.8)
    ax.set_xticks(xx); ax.set_xticklabels([f"{a}\n{b}" for a, b in allpairs], fontsize=10)
    ax.set_ylabel("delta dCor (Idasa - DMSO)")
    ax.set_title("(F) Specificity: only MDM2-CDKN1A gains dCor when loop opens;\n"
                 "others/controls drop (loss of cell-cycle covariation)", fontsize=12)

    cnt = {ln: (nd, ni) for ln, nd, ni, _ in ranking}
    nwt = cnt[wt_line]; nmut = cnt[mut_line]
    n_tot_d = int(data['DMSO']['cline'].size); n_tot_i = int(data['Idasanutlin']['cline'].size)
    fig.suptitle(f"Step 7 - Validation on REAL scRNA-seq (MIX-seq): CLOSED (DMSO) vs OPEN (Idasanutlin)\n"
                 f"[QC-passing singlets: DMSO {n_tot_d}, Idasa {n_tot_i} across 24 lines | "
                 f"WT {wt_line.split('_')[0]}: n={nwt[0]}/{nwt[1]}; "
                 f"MUT {mut_line.split('_')[0]}: n={nmut[0]}/{nmut[1]}  (DMSO/Idasa)]", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(str(FIG / "fig_step7_realdata.png"), dpi=120)
    print("Saved fig_step7_realdata.png")


if __name__ == "__main__":
    main()
