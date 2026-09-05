# -*- coding: utf-8 -*-
"""
p53-MDM2 model (MAIN) -- follows the Patterns 2021 template: a model is defined by
(propensities(x), a transitions/stoichiometry matrix, initial_state) and is simulated with
the authors' Gillespie engine ssa.py.

This is the classic p53-MDM2 oscillator formulation (Lev Bar-Or 2000; Geva-Zatorsky 2006;
Ciliberto/Tyson 2005; Proctor & Gray 2008): a delayed NEGATIVE feedback loop with
  - Michaelis-Menten (saturating, ultrasensitive) degradation of p53 by MDM2,
  - Hill (ultrasensitive) activation of MDM2 transcription by p53,
  - a time delay through the MDM2 mRNA intermediate -> enough to sustain oscillations.
The nonlinear propensities (Hill, MM) are valid hazard rates for the Markov-jump (Gillespie)
process, in the spirit of the "general propensity" used by Raharinirina et al. (2021).

Species (4) -- the scRNA-seq-observable mRNA pair is (p53_mRNA, Mdm2_mRNA):
    x[0] = p53_mRNA
    x[1] = p53        (protein)
    x[2] = Mdm2_mRNA
    x[3] = Mdm2       (protein)

Arm 1 (POSITIVE, transcription): p53 protein -> activates MDM2 mRNA transcription (Hill)  [R5]
Arm 2 (NEGATIVE, protein):       MDM2 protein -> degrades p53 (Michaelis-Menten)          [R4]
Nutlin-3: occupies the p53 pocket of MDM2 -> MDM2 cannot bind/degrade p53 -> CUT R4
          (k_degp53 -> k_degp53*(1-nutlin)). Result: p53 accumulates; p53 still drives MDM2
          transcription (R5 intact) so MDM2 mRNA rises, but cannot restrain p53 -> OPEN loop.
"""

import numpy as np

species = ('p53_mRNA', 'p53', 'Mdm2_mRNA', 'Mdm2')
IDX = {s: i for i, s in enumerate(species)}

# ---------------------------------------------------------------------------
# Parameters (units: /min, molecule counts). Tuned (see step 1) for sustained oscillations,
# a period of a few hours (matching experimental p53 pulses), with p53 & MDM2 on comparable scales.
# ---------------------------------------------------------------------------
DEFAULT_PARAMS = dict(
    # --- p53 mRNA and protein ---
    k_txn_p53      = 3.0,    # R0  p53 transcription (CONSTANT: p53 is regulated post-translationally)
    k_deg_p53mRNA  = 0.20,   # R1  p53 mRNA degradation
    k_tsl_p53      = 4.0,    # R2  p53 protein translation
    k_deg_p53      = 0.02,   # R3  basal p53 degradation (slow, MDM2-independent)
    # --- arm 2 (negative): MDM2-mediated p53 degradation, Michaelis-Menten ---  <-- NUTLIN blocks
    k_degp53       = 2.0,    # R4  maximal p53 degradation rate per MDM2
    Km_p53         = 120.0,  # R4  MM saturation constant in p53
    # --- arm 1 (positive): p53-activated MDM2 transcription, Hill ---
    k_txn_mdm2     = 6.0,    # R5  maximal MDM2 transcription rate
    Kd_mdm2        = 80.0,   # R5  p53 level for half-maximal activation (Hill)
    n_hill         = 4.0,    # R5  Hill coefficient (ultrasensitivity)
    k_txn_mdm2_0   = 0.05,   # R5  basal (leaky) MDM2 transcription
    # --- MDM2 mRNA and protein ---
    k_deg_mdm2mRNA = 0.10,   # R6  MDM2 mRNA degradation
    k_tsl_mdm2     = 1.2,    # R7  MDM2 protein translation
    k_deg_mdm2     = 0.30,   # R8  MDM2 protein degradation
)

# Initial (basal, unstressed) state: [p53_mRNA, p53, Mdm2_mRNA, Mdm2]
initial_state = np.array([15, 30, 20, 50])

# ---------------------------------------------------------------------------
# Stoichiometry matrix: transitions[:, j] = state change when reaction j fires.
# Species order: (p53_mRNA, p53, Mdm2_mRNA, Mdm2)
# ---------------------------------------------------------------------------
transitions = np.array([
    ( 1, 0, 0, 0),  # R0  0 -> p53_mRNA                 (p53 transcription)
    (-1, 0, 0, 0),  # R1  p53_mRNA -> 0                 (degradation)
    ( 0, 1, 0, 0),  # R2  p53_mRNA -> p53_mRNA + p53    (translation)
    ( 0,-1, 0, 0),  # R3  p53 -> 0                      (basal degradation)
    ( 0,-1, 0, 0),  # R4  p53 --(MDM2, MM)--> 0         (MDM2-mediated degradation) [ARM 2, Nutlin cuts]
    ( 0, 0, 1, 0),  # R5  0 --(p53, Hill)--> Mdm2_mRNA  (p53-activated transcription) [ARM 1]
    ( 0, 0,-1, 0),  # R6  Mdm2_mRNA -> 0                (degradation)
    ( 0, 0, 0, 1),  # R7  Mdm2_mRNA -> Mdm2_mRNA + Mdm2 (translation)
    ( 0, 0, 0,-1),  # R8  Mdm2 -> 0                     (protein degradation)
]).T

NUM_REACTIONS = transitions.shape[1]


def make_propensities(params=None, nutlin=0.0):
    """
    Return a propensities(x) function in the Patterns-2021 style.
    nutlin in [0,1]: Nutlin-3 efficacy; multiplies k_degp53 by (1-nutlin) (cuts the MDM2->p53 arm).
    """
    p = dict(DEFAULT_PARAMS)
    if params:
        p.update(params)
    kdeg_eff = p['k_degp53'] * (1.0 - nutlin)
    n = p['n_hill']

    def propensities(x):
        P = x[1]
        D = x[3]
        hill = P**n / (p['Kd_mdm2']**n + P**n)              # p53-activated MDM2 transcription
        mm   = D * P / (p['Km_p53'] + P)                    # MDM2-mediated p53 degradation (MM)
        return np.array([
            p['k_txn_p53'],                                 # R0
            p['k_deg_p53mRNA'] * x[0],                      # R1
            p['k_tsl_p53']     * x[0],                      # R2
            p['k_deg_p53']     * x[1],                      # R3
            kdeg_eff           * mm,                        # R4  [Nutlin: kdeg_eff]
            p['k_txn_mdm2'] * hill + p['k_txn_mdm2_0'],     # R5
            p['k_deg_mdm2mRNA']* x[2],                      # R6
            p['k_tsl_mdm2']    * x[2],                      # R7
            p['k_deg_mdm2']    * x[3],                      # R8
        ])
    return propensities


def ode_rhs(t, y, params=None, nutlin=0.0):
    """Mean-field (deterministic ODE) approximation of the same network -- for tuning/checks."""
    p = dict(DEFAULT_PARAMS)
    if params:
        p.update(params)
    kdeg_eff = p['k_degp53'] * (1.0 - nutlin)
    n = p['n_hill']
    mp, P, md, D = y
    hill = P**n / (p['Kd_mdm2']**n + P**n)
    mm   = D * P / (p['Km_p53'] + P)
    dmp = p['k_txn_p53'] - p['k_deg_p53mRNA'] * mp
    dP  = p['k_tsl_p53'] * mp - p['k_deg_p53'] * P - kdeg_eff * mm
    dmd = p['k_txn_mdm2'] * hill + p['k_txn_mdm2_0'] - p['k_deg_mdm2mRNA'] * md
    dD  = p['k_tsl_mdm2'] * md - p['k_deg_mdm2'] * D
    return [dmp, dP, dmd, dD]
