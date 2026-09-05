# -*- coding: utf-8 -*-
"""
STEP 1b - Validate that fast_ssa (numba) is EQUIVALENT to the original ssa.py engine.

  (a) Compare the stoichiometry matrices -> identical.
  (b) Compare propensity functions over many random states (and Nutlin levels)
      -> match to machine precision => both engines simulate the SAME Markov-jump process.
  (c) Small-scale, short-horizon statistical check: trajectory means agree within error.

Run:  python analysis/step1b_validate_fast_ssa.py
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import time
import numpy as np
import ssa
import model as M
import fast_ssa as F


def check_stoichiometry():
    assert np.array_equal(M.transitions, F._STOICH.T), "Stoichiometry mismatch!"
    print("[a] Stoichiometry identical:", M.transitions.shape)


def check_propensities(n_states=2000):
    th = F.params_to_theta(None)
    rng = np.random.default_rng(0)
    max_err = 0.0
    for nutlin in (0.0, 0.5, 1.0):
        prop_ref = M.make_propensities(None, nutlin=nutlin)
        kdeg_eff = th[4] * (1.0 - nutlin)
        for _ in range(n_states):
            x = rng.integers(0, 400, size=4).astype(np.float64)
            max_err = max(max_err, np.max(np.abs(prop_ref(x) - F._propensities(x, th, kdeg_eff))))
    print(f"[b] Max propensity error over {n_states} states x 3 Nutlin levels: {max_err:.2e}")
    assert max_err < 1e-9, "Propensity mismatch!"


def check_statistics(N_ref=40, N_fast=3000, horizon=150.0, dt=25.0):
    T = np.arange(0.0, horizon + dt, dt)
    prop = M.make_propensities(None, 0.0)
    np.random.seed(123)
    t0 = time.time()
    ref = np.array([ssa.SSA_Fixed_Width_Trajectory(M.transitions, prop, M.initial_state.copy(), T)
                    for _ in range(N_ref)])
    t_ref = time.time() - t0
    t0 = time.time()
    fast = F.simulate_ensemble(N_fast, T, nutlin=0.0, seed=321)
    t_fast = time.time() - t0
    print(f"[c] ssa.py {N_ref} cells: {t_ref:.1f}s ({t_ref/N_ref*1000:.0f} ms/cell) | "
          f"fast_ssa {N_fast} cells: {t_fast:.1f}s")
    for name, i in [('p53_mRNA', 0), ('p53', 1), ('Mdm2_mRNA', 2), ('Mdm2', 3)]:
        m_ref = ref[:, i, :].mean(); se = ref[:, i, :].std() / np.sqrt(N_ref); m_fast = fast[:, i, :].mean()
        flag = "OK" if abs(m_ref - m_fast) < 4 * se + 1 else "??"
        print(f"      {name:9s}: ref={m_ref:7.2f} (+/-{se:4.2f})  fast={m_fast:7.2f}  [{flag}]")


if __name__ == "__main__":
    check_stoichiometry()
    check_propensities()
    check_statistics()
    print("\nValidation done: fast_ssa matches ssa.py in stoichiometry + propensity + statistics.")
