# -*- coding: utf-8 -*-
"""Generate manuscript-quality standalone panels from Results 4.1-4.4.

Outputs are saved as both PNG (600 dpi) and PDF under figures/manuscript/.
The panel files intentionally do not contain panel letters.

Run:
    python analysis/manuscript_figure1_panels.py
"""
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
OUT = ROOT / "figures" / "manuscript"

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches

import model as M


COL = {
    "tp53_mrna": "#6b7280",
    "p53": "#d62728",
    "mdm2_mrna": "#1f77b4",
    "mdm2": "#2ca02c",
    "nutlin": "#b45309",
    "edge_pos": "#1f9d55",
    "edge_neg": "#c2410c",
    "closed": "#111827",
    "nutlin_line": "#b45309",
}


def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
    })


def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"saved {OUT / (stem + '.png')}")


def smooth(y, window=41):
    y = np.asarray(y, float)
    if len(y) < window:
        return y
    window = window if window % 2 else window + 1
    kernel = np.hanning(window)
    kernel /= kernel.sum()
    pad = window // 2
    return np.convolve(np.pad(y, pad, mode="edge"), kernel, mode="valid")


def normalize(y):
    y = np.asarray(y, float)
    lo, hi = np.nanmin(y), np.nanmax(y)
    return (y - lo) / (hi - lo + 1e-12)


def draw_gene_box(ax, xy, text, color, measured=False, width=1.55):
    x, y = xy
    box = patches.FancyBboxPatch(
        (x - width / 2, y - 0.28), width, 0.56,
        boxstyle="round,pad=0.04,rounding_size=0.08",
        linewidth=1.0, edgecolor=color, facecolor="white"
    )
    ax.add_patch(box)
    ax.text(x, y + 0.03, text, ha="center", va="center",
            fontsize=8.0, fontweight="bold", color=color, linespacing=0.9)
    if measured:
        ax.text(x, y - 0.43, "scRNA-seq observable", ha="center", va="center",
                fontsize=5.2, color="#374151")


def arrow(ax, start, end, color, lw=1.7, inhibit=False, **kwargs):
    if inhibit:
        ax.annotate("", xy=end, xytext=start,
                    arrowprops=dict(arrowstyle="-", color=color, lw=lw,
                                    shrinkA=8, shrinkB=13,
                                    connectionstyle=kwargs.get("connectionstyle", "arc3,rad=0")))
        x0, y0 = start
        x1, y1 = end
        vx, vy = x1 - x0, y1 - y0
        norm = (vx * vx + vy * vy) ** 0.5 + 1e-12
        px, py = -vy / norm, vx / norm
        xb, yb = x1 - 0.04 * vx / norm, y1 - 0.04 * vy / norm
        ax.plot([xb - 0.022 * px, xb + 0.022 * px],
                [yb - 0.022 * py, yb + 0.022 * py],
                color=color, lw=lw, solid_capstyle="round")
    else:
        ax.annotate("", xy=end, xytext=start,
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                    shrinkA=8, shrinkB=8, mutation_scale=12,
                                    connectionstyle=kwargs.get("connectionstyle", "arc3,rad=0")))


def panel_a():
    fig, ax = plt.subplots(figsize=(6.3, 2.9))
    ax.set_xlim(-0.15, 10.25)
    ax.set_ylim(0, 5)
    ax.axis("off")

    pos = {
        "tp53": (1.00, 3.65),
        "p53": (3.28, 3.65),
        "mdm2m": (6.28, 3.65),
        "mdm2": (8.92, 3.65),
    }
    draw_gene_box(ax, pos["tp53"], "TP53\nmRNA", COL["tp53_mrna"], measured=True, width=1.35)
    draw_gene_box(ax, pos["p53"], "p53\nprotein", COL["p53"], measured=False, width=1.35)
    draw_gene_box(ax, pos["mdm2m"], "MDM2\nmRNA", COL["mdm2_mrna"], measured=True, width=1.45)
    draw_gene_box(ax, pos["mdm2"], "MDM2\nprotein", COL["mdm2"], measured=False, width=1.45)

    arrow(ax, (1.73, 3.65), (2.58, 3.65), "#4b5563", lw=1.1)
    ax.text(2.15, 4.00, "translation", ha="center", fontsize=6.3, color="#4b5563")
    arrow(ax, (3.98, 3.65), (5.42, 3.65), COL["edge_pos"], lw=1.7)
    ax.text(4.72, 4.08, "transcriptional activation", ha="center",
            fontsize=6.5, color=COL["edge_pos"])
    arrow(ax, (7.08, 3.65), (8.08, 3.65), "#4b5563", lw=1.1)
    ax.text(7.58, 4.00, "translation", ha="center", fontsize=6.3, color="#4b5563")

    arrow(ax, (8.65, 3.22), (3.55, 3.15), COL["edge_neg"], lw=1.7, inhibit=True,
          connectionstyle="arc3,rad=-0.42")
    ax.text(6.10, 1.88, "post-translational inhibition: p53 degradation",
            ha="center", fontsize=6.8, color=COL["edge_neg"])

    nut = patches.FancyBboxPatch((4.95, 0.92), 2.40, 0.43,
                                 boxstyle="round,pad=0.04,rounding_size=0.08",
                                 linewidth=1.0, edgecolor=COL["nutlin"],
                                 facecolor="#fff7ed")
    ax.add_patch(nut)
    ax.text(6.15, 1.135, "Nutlin / idasanutlin", ha="center", va="center",
            fontsize=6.8, fontweight="bold", color=COL["nutlin"])
    ax.plot([5.91, 6.39], [2.30, 2.72], color=COL["nutlin"], lw=1.7)
    ax.plot([5.91, 6.39], [2.72, 2.30], color=COL["nutlin"], lw=1.7)
    arrow(ax, (6.15, 1.40), (6.15, 2.24), COL["nutlin"], lw=1.0)

    inset = patches.FancyBboxPatch((0.42, 0.42), 4.2, 0.95,
                                   boxstyle="round,pad=0.08,rounding_size=0.10",
                                   linewidth=0.9, edgecolor="#9ca3af",
                                   facecolor="#f9fafb")
    ax.add_patch(inset)
    ax.text(2.52, 1.14, "mRNA-only graph", ha="center", va="center",
            fontsize=7.0, fontweight="bold", color="#111827")
    ax.text(2.52, 0.84, "TP53 / p53 activity -> MDM2 mRNA", ha="center",
            fontsize=6.6, color=COL["edge_pos"])
    ax.text(2.52, 0.58, "no direct MDM2 mRNA -| TP53 mRNA edge", ha="center",
            fontsize=6.6, color=COL["edge_neg"])

    save(fig, "figure1A")


def panel_b(step1):
    tt = step1["tt"]
    y = step1["ode_0"]
    fig, ax = plt.subplots(figsize=(3.55, 2.65))
    ax.plot(tt, normalize(y[M.IDX["p53"]]), color=COL["p53"], lw=1.8, label="p53 protein")
    ax.plot(tt, normalize(y[M.IDX["Mdm2_mRNA"]]), color=COL["mdm2_mrna"], lw=1.8, label="MDM2 mRNA")
    ax.plot(tt, normalize(y[M.IDX["Mdm2"]]), color=COL["mdm2"], lw=1.8, label="MDM2 protein")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Normalized abundance")
    ax.set_title("Deterministic closed-loop dynamics")
    ax.set_xlim(0, 80)
    ax.legend(frameon=False, loc="upper right")
    ax.grid(True, color="#d1d5db", lw=0.5, alpha=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "figure1B")


def panel_c(step1):
    t = step1["T"]
    cells = step1["ssa_0"][:, M.IDX["p53"], :]
    fig, ax = plt.subplots(figsize=(3.55, 2.65))
    for i, trace in enumerate(cells):
        ax.plot(t, smooth(trace, 21), lw=1.1, alpha=0.9,
                color=plt.cm.tab10(i), label=f"cell {i + 1}")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("p53 protein count")
    ax.set_title("Stochastic single-cell p53 pulses")
    ax.legend(frameon=False, ncol=2, loc="upper right", handlelength=1.5)
    ax.grid(True, color="#d1d5db", lw=0.5, alpha=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "figure1C")


def panel_d(step1):
    tt = step1["tt"]
    closed = step1["ode_0"]
    nutlin = step1["ode_1"]
    fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.55), sharex=True)
    for ax, y, title in zip(axes, [closed, nutlin], ["Closed loop", "Nutlin / idasanutlin"]):
        p53_fc = np.log2((y[M.IDX["p53"]] + 1.0) / (y[M.IDX["p53"], 0] + 1.0))
        mdm2_fc = np.log2((y[M.IDX["Mdm2_mRNA"]] + 1.0) / (y[M.IDX["Mdm2_mRNA"], 0] + 1.0))
        ax.plot(tt, p53_fc, color=COL["p53"], lw=1.9, label="p53 protein")
        ax.plot(tt, mdm2_fc, color=COL["mdm2_mrna"], lw=1.9, label="MDM2 mRNA")
        ax.set_title(title)
        ax.set_xlabel("Time (min)")
        ax.set_xlim(0, 300)
        ax.grid(True, color="#d1d5db", lw=0.5, alpha=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("log2 fold-change from t = 0")
    axes[1].legend(frameon=False, loc="lower right")
    save(fig, "figure1D")


def panel_e(step2):
    dose = step2["dose"]
    fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.25), sharex=True)
    specs = [
        ("p53_mean", "Mean p53 protein", COL["p53"], "o"),
        ("mdm2_mrna_mean", "Mean MDM2 mRNA", COL["mdm2_mrna"], "s"),
        ("p53_osc", "Oscillation index", COL["mdm2"], "^"),
    ]
    for ax, (key, ylabel, color, marker) in zip(axes, specs):
        ax.plot(dose, step2[key], marker=marker, ms=4.2, lw=1.8, color=color)
        ax.set_xlabel("Nutlin efficacy")
        ax.set_ylabel(ylabel)
        ax.set_xlim(-0.03, 1.03)
        ax.grid(True, color="#d1d5db", lw=0.5, alpha=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    save(fig, "figure1E")


def main():
    style()
    step1 = np.load(DATA / "cache_step1.npz", allow_pickle=True)
    step2 = np.load(DATA / "cache_step2.npz", allow_pickle=True)
    panel_a()
    panel_b(step1)
    panel_c(step1)
    panel_d(step1)
    panel_e(step2)


if __name__ == "__main__":
    main()
