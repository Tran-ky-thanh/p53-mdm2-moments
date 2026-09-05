# -*- coding: utf-8 -*-
"""
STEP 1 - Single-cell dynamics of the p53-MDM2 loop, using the ORIGINAL Gillespie
engine (ssa.py) from Raharinirina et al. (Patterns 2021).

Deterministic ODE (damped) vs stochastic SSA (noise-sustained oscillations), closed loop
vs NUTLIN-3. Simulation results are cached to data/cache_step1.npz; re-plotting only reloads.

Run:            python analysis/step1_single_cell_dynamics.py
Force recompute: RECOMPUTE=1 python analysis/step1_single_cell_dynamics.py
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
FIG = ROOT / "figures"; DATA = ROOT / "data"

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import ssa
import model as M
import plotstyle
from cache import load_or_compute

HORIZON = 1500.0
N_CELLS = 4
CONDS = [(0.0, "CLOSED loop (no drug)"), (1.0, "NUTLIN-3 (MDM2->p53 arm cut)")]


def compute():
    out = {}
    tt = np.linspace(0, HORIZON, 3000)
    T = np.arange(0.0, HORIZON, 2.0)
    out['tt'] = tt; out['T'] = T
    for col, (nutlin, _) in enumerate(CONDS):
        sol = solve_ivp(lambda t, y: M.ode_rhs(t, y, None, nutlin), [0, HORIZON],
                        M.initial_state.astype(float), max_step=1.0, rtol=1e-7, atol=1e-9,
                        dense_output=True)
        out[f'ode_{col}'] = sol.sol(tt)
        np.random.seed(10 + col)
        prop = M.make_propensities(None, nutlin=nutlin)
        cells = [ssa.SSA_Fixed_Width_Trajectory(M.transitions, prop, M.initial_state.copy(), T)
                 for _ in range(N_CELLS)]
        out[f'ssa_{col}'] = np.array(cells)
    return out


def plot(D):
    plotstyle.apply()
    fig, axs = plt.subplots(2, 2, figsize=(14, 8))
    tt = D['tt']; T = D['T']
    for col, (nutlin, title) in enumerate(CONDS):
        Y = D[f'ode_{col}']
        ax = axs[0, col]
        ax.plot(tt, Y[1], label='p53 protein', color='C0')
        ax.plot(tt, Y[2], label='MDM2 mRNA', color='C1')
        ax.plot(tt, Y[3], label='MDM2 protein', color='C2', alpha=0.6)
        ax.set_title(f"ODE (deterministic) - {title}", fontsize=10)
        ax.set_ylabel('molecule count'); ax.legend(fontsize=8)
        ax = axs[1, col]
        for c in D[f'ssa_{col}']:
            ax.plot(T, c[1], color='C0', lw=0.7, alpha=0.8)
            ax.plot(T, c[2], color='C1', lw=0.7, alpha=0.8)
        ax.plot([], [], color='C0', label='p53 protein'); ax.plot([], [], color='C1', label='MDM2 mRNA')
        ax.set_title(f"SSA, {N_CELLS} cells (ssa.py engine) - {title}", fontsize=10)
        ax.set_xlabel('time (min)'); ax.set_ylabel('molecule count'); ax.legend(fontsize=8)
    fig.suptitle(f"Step 1 - p53-MDM2 negative feedback: oscillations (closed loop) vs OPEN loop (Nutlin-3)  "
                 f"[deterministic ODE + {N_CELLS} SSA cells per condition, horizon {HORIZON:.0f} min]",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(str(FIG / "fig_step1_single_cell.png"), dpi=120)
    print("Saved fig_step1_single_cell.png")


def main():
    D = load_or_compute(DATA / "cache_step1.npz", compute)
    plot(D)


if __name__ == "__main__":
    main()
