# -*- coding: utf-8 -*-
"""STEP 9b - Synthetic p53-MDM2-CDKN1A data for Figure 13D.

The six-species Gillespie model contains mRNA and protein for all three genes.
Nutlin efficacy is fixed at 1, so the MDM2-mediated p53 degradation arm is cut,
while p53 still activates transcription of both MDM2 and CDKN1A.

The complete simulated observation array is saved for later analyses:
    data/cache_step9_three_gene_nutlin.npz

Run:            python analysis/step9b_three_gene_nutlin.py
Force recompute: RECOMPUTE=1 python analysis/step9b_three_gene_nutlin.py
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"

import numpy as np
import fast_ssa_three_gene as F3
import model_three_gene as M3
import directional as D
from cache import load_or_compute

N_CELLS = 4000
TIMES = np.arange(0.0, 81.0, 4.0)
NUTLIN = 1.0
SEED = 93
CACHE = DATA / "cache_step9_three_gene_nutlin.npz"


def compute():
    trajectories = F3.simulate_ensemble(N_CELLS, TIMES, nutlin=NUTLIN, seed=SEED)
    md = trajectories[:, M3.IDX['Mdm2_mRNA'], :]
    cd = trajectories[:, M3.IDX['CDKN1A_mRNA'], :]
    mp = trajectories[:, M3.IDX['p53_mRNA'], :]
    pp = trajectories[:, M3.IDX['p53'], :]

    # At each time, condition either on the observed TP53 transcript or on recent
    # p53-protein history. The latter is the actual common regulator in the model.
    stats = np.zeros((len(TIMES), 3))
    for i in range(1, len(TIMES)):
        lo = max(0, i - 4)
        stats[i] = [np.corrcoef(md[:, i], cd[:, i])[0, 1],
                    D.partial_corr(md[:, i], cd[:, i], mp[:, i]),
                    D.partial_corr(md[:, i], cd[:, i], pp[:, lo:i + 1])]
    return dict(times=TIMES, trajectories=trajectories, stats=stats,
                species=np.array(M3.species), nutlin=np.array([NUTLIN]),
                n_cells=np.array([N_CELLS]), seed=np.array([SEED]))


def main():
    out = load_or_compute(CACHE, compute)
    print(f"Saved simulation: {CACHE}")
    print("species:", list(out['species']))
    for i, t in enumerate(out['times']):
        if t > 0:
            print(f"t={t:4.0f}: marginal={out['stats'][i,0]:+.3f}, "
                  f"partial|TP53mRNA={out['stats'][i,1]:+.3f}, "
                  f"partial|p53protein-history={out['stats'][i,2]:+.3f}")


if __name__ == "__main__":
    main()
