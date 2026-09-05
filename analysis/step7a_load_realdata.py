# -*- coding: utf-8 -*-
"""
STEP 7a - Load REAL scRNA-seq data (MIX-seq, McFarland et al. Nat Commun 2020).

Closed vs open loop on real data:
  - DMSO_6hr_expt1        = vehicle control  -> "CLOSED loop" (intact p53-MDM2)
  - Idasanutlin_6hr_expt1 = MDM2 inhibitor (Nutlin family, RG7388) -> "OPEN loop"
Pool expt1 = 24 cell lines; Nutlin activates p53 only in TP53 wild-type (WT) lines.

This loader reads matrix.mtx (genes x cells), keeps QC-passing 'normal' singlets, extracts a
gene PANEL (p53 targets + control genes), per-cell total counts (for library-size normalization),
and cell-line identity; saves a compact .npz for step 7b.

Set DATA_DIR to the folder that contains the MIX-seq experiment subfolders.
Run:  python analysis/step7a_load_realdata.py
"""
import sys, pathlib, os, csv
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

import numpy as np
import scipy.io as sio

# Folder holding the MIX-seq experiment subfolders (edit if needed):
DATA_DIR = str((ROOT.parent / "10298696").resolve())
# Idasanutlin (2.5 uM) was assayed at two timepoints in expt1: 6 h and 24 h.
TIMEPOINTS = ["6hr", "24hr"]

P53_TARGETS = ["MDM2", "CDKN1A", "FDXR", "TP53I3", "BTG2", "GDF15", "SESN1",
               "RRM2B", "DDB2", "BAX", "ZMAT3", "PHLDA3", "TNFRSF10B", "AEN", "TP53"]
CONTROL_GENES = ["ACTB", "GAPDH", "B2M", "TUBB", "RPL13A", "PGK1", "HPRT1",
                 "PPIA", "TBP", "GUSB", "VIM", "EEF1A1"]
PANEL = P53_TARGETS + CONTROL_GENES


def load_gene_index(genes_tsv):
    sym2rows = {}
    with open(genes_tsv) as f:
        for i, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            sym = parts[1] if len(parts) > 1 else parts[0]
            sym2rows.setdefault(sym, []).append(i)
    return sym2rows


def load_condition(folder):
    path = os.path.join(DATA_DIR, folder)
    print(f"  reading {folder} ...", flush=True)
    Mx = sio.mmread(os.path.join(path, "matrix.mtx")).tocsr()   # genes x cells
    n_genes, n_cells = Mx.shape
    tot = np.asarray(Mx.sum(axis=0)).ravel()
    sym2rows = load_gene_index(os.path.join(path, "genes.tsv"))
    with open(os.path.join(path, "barcodes.tsv")) as f:
        barcodes = [l.strip() for l in f]
    q, cl = {}, {}
    with open(os.path.join(path, "classifications.csv")) as f:
        for row in csv.DictReader(f):
            q[row["barcode"]] = row["cell_quality"]; cl[row["barcode"]] = row["singlet_ID"]
    quality = np.array([q.get(bc, "NA") for bc in barcodes])
    cline = np.array([cl.get(bc, "NA") for bc in barcodes])
    panel_counts = {}
    for sym in PANEL:
        rows = sym2rows.get(sym, [])
        if not rows:
            continue
        r = rows[0] if len(rows) == 1 else rows[int(np.argmax([Mx[rr].sum() for rr in rows]))]
        panel_counts[sym] = np.asarray(Mx[r].todense()).ravel().astype(np.float64)
    return dict(counts=panel_counts, tot=tot, quality=quality, cline=cline, n_cells=n_cells)


def build_panel(tp):
    """Build the (DMSO, Idasanutlin) panel for one timepoint tp ('6hr' or '24hr')."""
    out = {}
    for cond in ("DMSO", "Idasanutlin"):
        d = load_condition(f"{cond}_{tp}_expt1")
        keep = d["quality"] == "normal"
        print(f"  [{tp}] {cond}: {d['n_cells']} cells -> {keep.sum()} 'normal' singlets")
        genes_found = [g for g in PANEL if g in d["counts"]]
        out[f"{cond}__X"] = np.stack([d["counts"][g][keep] for g in genes_found], axis=1)
        out[f"{cond}__genes"] = np.array(genes_found)
        out[f"{cond}__tot"] = d["tot"][keep]
        out[f"{cond}__cline"] = d["cline"][keep]
    return out


def main():
    for tp in TIMEPOINTS:
        out = build_panel(tp)
        # 6 h stays as the canonical panel; 24 h gets a suffixed file
        fname = "data_step7_panel.npz" if tp == "6hr" else f"data_step7_panel_{tp}.npz"
        np.savez_compressed(str(DATA / fname), **out)
        print(f"Saved data/{fname}")
    print("Genes not found:", [g for g in PANEL if g not in set(out["DMSO__genes"])] or "(none)")


if __name__ == "__main__":
    main()
