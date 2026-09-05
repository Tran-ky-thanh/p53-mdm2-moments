# -*- coding: utf-8 -*-
"""
STEP 2 - Nutlin-3 DOSE response (p53 up, MDM2 mRNA up, oscillation down).
Simulation cached to data/cache_step2.npz; re-plotting only reloads.

Run:            python analysis/step2_nutlin_dose_response.py
Force recompute: RECOMPUTE=1 python analysis/step2_nutlin_dose_response.py
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
FIG = ROOT / "figures"; DATA = ROOT / "data"

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fast_ssa as F
import model as M
import plotstyle
from cache import load_or_compute

N_CELLS = 1500
HORIZON = 1600.0
DT = 10.0
DOSES = [0.0, 0.25, 0.5, 0.75, 0.9, 1.0]
LATE_FROM = 600.0


def oscillation_index(traj_late):
    m = traj_late.mean(axis=1); s = traj_late.std(axis=1)
    return (s / np.maximum(m, 1e-9)).mean()


def compute():
    T = np.arange(0.0, HORIZON, DT)
    late = T >= LATE_FROM
    dose, p53m, mdm2m, osc = [], [], [], []
    for nut in DOSES:
        data = F.simulate_ensemble(N_CELLS, T, nutlin=nut, seed=1000 + int(nut * 100))
        p53 = data[:, M.IDX['p53'], :][:, late]
        mdm = data[:, M.IDX['Mdm2_mRNA'], :][:, late]
        dose.append(nut); p53m.append(p53.mean()); mdm2m.append(mdm.mean()); osc.append(oscillation_index(p53))
        print(f"nutlin={nut:.2f}: p53 mean={p53.mean():7.1f}  MDM2 mRNA mean={mdm.mean():6.1f}  "
              f"osc(CV)={osc[-1]:.3f}")
    return dict(dose=np.array(dose), p53_mean=np.array(p53m),
                mdm2_mrna_mean=np.array(mdm2m), p53_osc=np.array(osc))


def plot(D):
    plotstyle.apply()
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
    axs[0].plot(D['dose'], D['p53_mean'], 'o-', color='C0')
    axs[0].set_title("p53 protein accumulates with Nutlin dose"); axs[0].set_ylabel("mean p53 (molecules)")
    axs[1].plot(D['dose'], D['mdm2_mrna_mean'], 's-', color='C1')
    axs[1].set_title("MDM2 mRNA rises (transcription arm intact)"); axs[1].set_ylabel("mean MDM2 mRNA")
    axs[2].plot(D['dose'], D['p53_osc'], '^-', color='C2')
    axs[2].set_title("p53 oscillation dies as the loop opens"); axs[2].set_ylabel("oscillation index (temporal CV)")
    for ax in axs:
        ax.set_xlabel("Nutlin-3 efficacy  (0 = no drug, 1 = full block)"); ax.grid(alpha=0.3)
    fig.suptitle(f"Step 2 - Nutlin-3 dose response ({N_CELLS}-cell populations, window t>={LATE_FROM:.0f} min)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(str(FIG / "fig_step2_dose_response.png"), dpi=120)
    print("Saved fig_step2_dose_response.png")


def main():
    D = load_or_compute(DATA / "cache_step2.npz", compute)
    plot(D)


if __name__ == "__main__":
    main()
