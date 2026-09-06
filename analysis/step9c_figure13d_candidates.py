# -*- coding: utf-8 -*-
"""STEP 9c - Draw four candidate replacements for Figure 13D.

All panels are derived from data/cache_step9_three_gene_nutlin.npz. This step
does not rerun the six-species Gillespie simulation. It saves four standalone
PNGs and a self-contained HTML comparison page.

Run: python analysis/step9c_figure13d_candidates.py
"""
import base64
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
OUT = ROOT / "figures" / "figure13D_candidates"
REPORT = ROOT / "report" / "figure13D_candidates.html"

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import plotstyle

CACHE = DATA / "cache_step9_three_gene_nutlin.npz"


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved", path)
    return path


def option1_scatter(z, idx):
    md = z["snapshot_mdm2_mRNA"]
    cd = z["snapshot_cdkn1a_mRNA"]
    p53 = z["snapshot_p53_exposure"][:, 0]
    r = np.corrcoef(md, cd)[0, 1]
    fig, ax = plt.subplots(figsize=(8.6, 7.2))
    pts = ax.scatter(md, cd, c=p53, cmap="viridis", s=16, alpha=.45,
                     edgecolors="none", rasterized=True)
    slope, intercept = np.polyfit(md, cd, 1)
    xx = np.linspace(md.min(), np.percentile(md, 99.5), 100)
    ax.plot(xx, intercept + slope * xx, color="crimson", lw=2.5,
            label=f"marginal trend (r = {r:+.3f})")
    cb = fig.colorbar(pts, ax=ax)
    cb.set_label("integrated p53-protein transcriptional input")
    ax.set_xlabel("MDM2 mRNA (molecules)")
    ax.set_ylabel("CDKN1A mRNA (molecules)")
    ax.set_title("Option 1 — Cell-level target co-expression under Nutlin\n"
                 "both targets rise with accumulated p53 activity")
    ax.legend(loc="upper left")
    ax.grid(alpha=.25)
    return save(fig, "fig13D_option1_scatter.png")


def option2_trajectories(z, idx):
    x = z["trajectories"]
    t = z["times"]
    curves = {
        "p53 protein": x[:, idx["p53"], :].mean(axis=0),
        "MDM2 mRNA": x[:, idx["Mdm2_mRNA"], :].mean(axis=0),
        "CDKN1A mRNA": x[:, idx["CDKN1A_mRNA"], :].mean(axis=0),
    }
    fig, ax = plt.subplots(figsize=(8.6, 7.2))
    colors = ["C3", "C0", "C2"]
    for (label, y), color in zip(curves.items(), colors):
        logfc = np.log2((y + 1) / (y[0] + 1))
        ax.plot(t, logfc, "o-", lw=2.6, ms=5, color=color, label=label)
    ax.axhline(0, color="gray", lw=.8)
    ax.set_xlabel("time after Nutlin (min)")
    ax.set_ylabel("log2 fold-change from t = 0")
    ax.set_title("Option 2 — Open-loop response dynamics\n"
                 "p53 accumulates; MDM2 and CDKN1A transcription remain active")
    ax.legend()
    ax.grid(alpha=.3)
    return save(fig, "fig13D_option2_trajectories.png")


def option3_heatmap(z, idx):
    exposure = z["snapshot_p53_exposure"][:, 0]
    values = np.c_[exposure, z["snapshot_mdm2_mRNA"], z["snapshot_cdkn1a_mRNA"]]
    order = np.argsort(exposure)
    take = order[np.linspace(0, len(order) - 1, 700).astype(int)]
    V = values[take].astype(float)
    lo = np.percentile(V, 1, axis=0); hi = np.percentile(V, 99, axis=0)
    V = np.clip((V - lo) / np.maximum(hi - lo, 1e-12), 0, 1)
    fig, ax = plt.subplots(figsize=(8.0, 7.2))
    im = ax.imshow(V, aspect="auto", cmap="magma", interpolation="nearest")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["p53 activity", "MDM2 mRNA", "CDKN1A mRNA"])
    ax.set_ylabel("cells sorted by increasing p53 activity")
    ax.set_yticks([0, len(V)-1]); ax.set_yticklabels(["low", "high"])
    ax.set_title("Option 3 — Cell-state heterogeneity under Nutlin\n"
                 "higher p53 exposure co-activates both target genes")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("within-feature scaled value")
    return save(fig, "fig13D_option3_heatmap.png")


def option4_response_path(z, idx):
    x = z["trajectories"]
    t = z["times"]
    md = x[:, idx["Mdm2_mRNA"], :].mean(axis=0)
    cd = x[:, idx["CDKN1A_mRNA"], :].mean(axis=0)
    points = np.c_[md, cd]
    segments = np.stack([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap="plasma", norm=plt.Normalize(t.min(), t.max()))
    lc.set_array(t[:-1]); lc.set_linewidth(3)
    fig, ax = plt.subplots(figsize=(8.6, 7.2))
    ax.add_collection(lc)
    sc = ax.scatter(md, cd, c=t, cmap="plasma", s=55, zorder=3)
    offsets = {0: (6, 5), 20: (6, 5), 40: (-42, 8), 60: (-42, -17), 80: (7, 9)}
    for target in [0, 20, 40, 60, 80]:
        k = int(np.argmin(np.abs(t - target)))
        ax.annotate(f"{t[k]:.0f} min", (md[k], cd[k]), xytext=offsets[target],
                    textcoords="offset points", fontsize=9)
    dx = max(md.max()-md.min(), 1e-6); dy = max(cd.max()-cd.min(), 1e-6)
    ax.set_xlim(md.min() - .05*dx, md.max() + .08*dx)
    ax.set_ylim(cd.min() - .05*dy, cd.max() + .08*dy)
    ax.set_xlabel("population mean MDM2 mRNA")
    ax.set_ylabel("population mean CDKN1A mRNA")
    ax.set_title("Option 4 — Population response path under Nutlin\n"
                 "p53 drives the two targets along a shared expression trajectory")
    cb = fig.colorbar(sc, ax=ax); cb.set_label("time after Nutlin (min)")
    ax.grid(alpha=.3)
    return save(fig, "fig13D_option4_response_path.png")


def img64(path):
    return base64.b64encode(path.read_bytes()).decode("ascii")


def write_html(paths):
    descriptions = [
        ("1. Cell-level scatter — recommended",
         "Each point is one simulated cell sampled once between 4 and 80 min. Colour is the integrated p53-protein input. The diagonal cloud shows that cells with stronger accumulated p53 activity express both MDM2 and CDKN1A more strongly. This connects directly to panel E: conditioning on p53 removes this common-driver correlation."),
        ("2. Mean response trajectories",
         "p53 protein rises after Nutlin because MDM2 can no longer degrade it. MDM2 mRNA and CDKN1A mRNA still rise because the p53-to-target transcriptional arms remain intact. This is the clearest picture of the open-loop biology, but it explains partial correlation less directly."),
        ("3. p53-sorted cell heatmap",
         "Cells are ordered by accumulated p53 activity. The coordinated transition from low to high MDM2 and CDKN1A shows treatment-response heterogeneity and common p53 control. It is compact and intuitive, although exact effect sizes are harder to read."),
        ("4. MDM2-CDKN1A response path",
         "The population mean moves through MDM2-CDKN1A expression space after Nutlin. The diagonal path shows their coordinated transcriptional response to p53 and makes the temporal order visible. It summarizes the population and hides individual-cell heterogeneity."),
    ]
    cards = []
    for path, (title, body) in zip(paths, descriptions):
        cards.append(f'<section><h2>{title}</h2><img src="data:image/png;base64,{img64(path)}"><p>{body}</p></section>')
    html = """<!doctype html><html><head><meta charset="utf-8"><title>Figure 13D candidates</title>
<style>body{font:18px/1.55 Arial,sans-serif;max-width:1500px;margin:30px auto;color:#1f2937}h1{text-align:center}main{display:grid;grid-template-columns:1fr 1fr;gap:28px}section{border:1px solid #d1d5db;border-radius:10px;padding:18px;background:#fff;box-shadow:0 2px 8px #0001}img{width:100%;height:auto}h2{font-size:22px;margin:0 0 8px}p{margin:10px 5px}@media(max-width:900px){main{grid-template-columns:1fr}}</style></head><body>
<h1>Candidate replacements for Figure 13D</h1><p>All candidates use the saved 4000-cell, six-species p53-MDM2-CDKN1A simulation with Nutlin efficacy = 1. No simulation was rerun.</p><main>""" + "".join(cards) + "</main></body></html>"
    REPORT.write_text(html, encoding="utf-8")
    print("Wrote", REPORT)


def main():
    plotstyle.apply()
    with np.load(CACHE, allow_pickle=True) as raw:
        z = {k: raw[k] for k in raw.files}
    idx = {str(name): i for i, name in enumerate(z["species"])}
    paths = [option1_scatter(z, idx), option2_trajectories(z, idx),
             option3_heatmap(z, idx), option4_response_path(z, idx)]
    write_html(paths)


if __name__ == "__main__":
    main()
