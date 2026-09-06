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
SNAPSHOT_SEED = 101
CACHE = DATA / "cache_step9_three_gene_nutlin.npz"


def asynchronous_snapshot(trajectories, times):
    """Take one observation per cell across heterogeneous response phases.

    The conditioning variables are two exponentially weighted summaries of the
    simulated p53-protein trajectory. They use the exact Hill input and mRNA
    decay rate of each target, so together they represent the shared p53 driver
    seen by MDM2 and CDKN1A before the sampled snapshot.
    """
    rng = np.random.default_rng(SNAPSHOT_SEED)
    n = trajectories.shape[0]
    obs_idx = rng.integers(1, len(times), n)
    cells = np.arange(n)
    md = trajectories[cells, M3.IDX['Mdm2_mRNA'], obs_idx]
    cd = trajectories[cells, M3.IDX['CDKN1A_mRNA'], obs_idx]
    protein = trajectories[:, M3.IDX['p53'], :]
    p = M3.DEFAULT_PARAMS
    hm = protein**p['n_hill'] / (p['Kd_mdm2']**p['n_hill'] + protein**p['n_hill'])
    hc = protein**p['n_cdkn1a'] / (p['Kd_cdkn1a']**p['n_cdkn1a'] + protein**p['n_cdkn1a'])
    dt = float(times[1] - times[0])
    exposure = np.zeros((n, 2))
    for cell, stop in enumerate(obs_idx):
        used = np.arange(stop + 1)
        age = times[stop] - times[used]
        exposure[cell, 0] = np.sum(hm[cell, used] * np.exp(-p['k_deg_mdm2mRNA'] * age)) * dt
        exposure[cell, 1] = np.sum(hc[cell, used] * np.exp(-p['k_deg_cdkn1a_mRNA'] * age)) * dt
    summary = np.array([np.corrcoef(md, cd)[0, 1], D.partial_corr(md, cd, exposure)])
    return dict(snapshot_time_index=obs_idx, snapshot_time=times[obs_idx],
                snapshot_mdm2_mRNA=md, snapshot_cdkn1a_mRNA=cd,
                snapshot_p53_exposure=exposure, snapshot_stats=summary)


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
    out = dict(times=TIMES, trajectories=trajectories, stats=stats,
               species=np.array(M3.species), nutlin=np.array([NUTLIN]),
               n_cells=np.array([N_CELLS]), seed=np.array([SEED]),
               snapshot_seed=np.array([SNAPSHOT_SEED]))
    out.update(asynchronous_snapshot(trajectories, TIMES))
    return out


def main():
    out = load_or_compute(CACHE, compute)
    if 'snapshot_stats' not in out:
        # Upgrade an existing cache without rerunning the expensive SSA simulation.
        out.update(asynchronous_snapshot(out['trajectories'], out['times']))
        out['snapshot_seed'] = np.array([SNAPSHOT_SEED])
        np.savez_compressed(CACHE, **out)
        print(f"[cache] added asynchronous snapshot to {CACHE.name}")
    print(f"Saved simulation: {CACHE}")
    print("species:", list(out['species']))
    for i, t in enumerate(out['times']):
        if t > 0:
            print(f"t={t:4.0f}: marginal={out['stats'][i,0]:+.3f}, "
                  f"partial|TP53mRNA={out['stats'][i,1]:+.3f}, "
                  f"partial|p53protein-history={out['stats'][i,2]:+.3f}")
    print(f"asynchronous snapshot: marginal={out['snapshot_stats'][0]:+.3f}, "
          f"partial(MDM2,CDKN1A|p53)={out['snapshot_stats'][1]:+.3f}")


if __name__ == "__main__":
    main()
