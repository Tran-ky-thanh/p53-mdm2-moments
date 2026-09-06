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
import model_three_gene as M3
import directional as D
import association as A
from cache import load_or_compute

N_CELLS = 4000
TIMES = np.arange(0.0, 81.0, 4.0)
NUTLIN = 1.0
SEED = 93
SNAPSHOT_SEED = 101
SPLATTER_SEED = 202
SPLATTER_LIB_SCALE = 0.35
SPLATTER_BCV = 0.40
SPLATTER_DROPOUT_MID = 1.0
SPLATTER_DROPOUT_SHAPE = -1.0
CACHE = DATA / "cache_step9_three_gene_nutlin.npz"


def simulation_metadata():
    """Self-describing metadata stored beside the trajectories."""
    names = np.array(list(M3.DEFAULT_PARAMS), dtype=str)
    values = np.array([M3.DEFAULT_PARAMS[name] for name in names], dtype=float)
    return dict(initial_state=M3.initial_state.astype(np.int64),
                parameter_names=names, parameter_values=values,
                model_version=np.array([2], dtype=np.int64))


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


def splatter_snapshot(snapshot_mdm2, snapshot_cdkn1a, exposure):
    """Apply Splatter-style technical noise to the synthetic snapshot targets."""
    rng = np.random.default_rng(SPLATTER_SEED)
    clean_counts = np.column_stack([snapshot_mdm2, snapshot_cdkn1a])
    noisy_counts = A.splatter_noise(clean_counts, rng, lib_scale=SPLATTER_LIB_SCALE,
                                    bcv=SPLATTER_BCV,
                                    dropout_mid=SPLATTER_DROPOUT_MID,
                                    dropout_shape=SPLATTER_DROPOUT_SHAPE)
    noisy_mdm2 = noisy_counts[:, 0]
    noisy_cdkn1a = noisy_counts[:, 1]
    stats = np.array([np.corrcoef(noisy_mdm2, noisy_cdkn1a)[0, 1],
                      D.partial_corr(noisy_mdm2, noisy_cdkn1a, exposure)])
    zero_rates = np.array([(noisy_mdm2 == 0).mean(),
                           (noisy_cdkn1a == 0).mean(),
                           (noisy_counts == 0).any(axis=1).mean(),
                           (noisy_counts == 0).all(axis=1).mean()], dtype=float)
    return dict(snapshot_splatter_counts=noisy_counts,
                snapshot_splatter_stats=stats,
                snapshot_splatter_zero_rates=zero_rates,
                snapshot_splatter_zero_rate_names=np.array([
                    "MDM2_zero_fraction",
                    "CDKN1A_zero_fraction",
                    "any_target_zero_fraction",
                    "both_targets_zero_fraction"
                ], dtype=str),
                snapshot_splatter_seed=np.array([SPLATTER_SEED], dtype=np.int64),
                snapshot_splatter_params=np.array([
                    SPLATTER_LIB_SCALE,
                    SPLATTER_BCV,
                    SPLATTER_DROPOUT_MID,
                    SPLATTER_DROPOUT_SHAPE,
                ], dtype=float),
                snapshot_splatter_param_names=np.array([
                    "lib_scale", "bcv", "dropout_mid", "dropout_shape"
                ], dtype=str))


def compute():
    # Import the numba engine only when trajectories actually need recomputing;
    # cached plotting and metadata upgrades therefore need only NumPy.
    import fast_ssa_three_gene as F3
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
    out.update(splatter_snapshot(out['snapshot_mdm2_mRNA'],
                                 out['snapshot_cdkn1a_mRNA'],
                                 out['snapshot_p53_exposure']))
    out.update(simulation_metadata())
    return out


def main():
    out = load_or_compute(CACHE, compute)
    upgraded = False
    if 'snapshot_stats' not in out:
        # Upgrade an existing cache without rerunning the expensive SSA simulation.
        out.update(asynchronous_snapshot(out['trajectories'], out['times']))
        out['snapshot_seed'] = np.array([SNAPSHOT_SEED])
        upgraded = True
    if 'snapshot_splatter_stats' not in out or 'snapshot_splatter_zero_rates' not in out:
        out.update(splatter_snapshot(out['snapshot_mdm2_mRNA'],
                                     out['snapshot_cdkn1a_mRNA'],
                                     out['snapshot_p53_exposure']))
        upgraded = True
    for key, value in simulation_metadata().items():
        if key not in out:
            out[key] = value
            upgraded = True
    if upgraded:
        np.savez_compressed(CACHE, **out)
        print(f"[cache] added analysis metadata to {CACHE.name}; trajectories were not recomputed")
    print(f"Saved simulation: {CACHE}")
    print("species:", list(out['species']))
    for i, t in enumerate(out['times']):
        if t > 0:
            print(f"t={t:4.0f}: marginal={out['stats'][i,0]:+.3f}, "
                  f"partial|TP53mRNA={out['stats'][i,1]:+.3f}, "
                  f"partial|p53protein-history={out['stats'][i,2]:+.3f}")
    print(f"asynchronous snapshot: marginal={out['snapshot_stats'][0]:+.3f}, "
          f"partial(MDM2,CDKN1A|p53)={out['snapshot_stats'][1]:+.3f}")
    print(f"splatter snapshot: marginal={out['snapshot_splatter_stats'][0]:+.3f}, "
          f"partial(MDM2,CDKN1A|p53)={out['snapshot_splatter_stats'][1]:+.3f}")
    print("splatter zero fractions:",
          {str(k): round(float(v), 4) for k, v in zip(out['snapshot_splatter_zero_rate_names'],
                                                      out['snapshot_splatter_zero_rates'])})


if __name__ == "__main__":
    main()
