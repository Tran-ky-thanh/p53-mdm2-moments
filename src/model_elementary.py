# -*- coding: utf-8 -*-
"""
ALTERNATIVE p53-MDM2 model (elementary mass-action) -- same Patterns 2021 template:
a model is (propensities(x), transitions/stoichiometry matrix, initial_state).

Reaction structure based on the stochastic model of Proctor & Gray 2008 (BMC Syst Biol 2:75):
a p53-MDM2 NEGATIVE feedback loop, with an MDM2 mRNA delay and p53 degradation via an Mdm2
complex. This variant is kept for reference; the main model.py uses MM/Hill propensities
instead, which give a more balanced fixed point (see the report / step 1).

Species (5):
    x[0] = p53_mRNA
    x[1] = p53        (free protein)
    x[2] = Mdm2_mRNA
    x[3] = Mdm2       (free protein)
    x[4] = Mdm2_p53   (complex)

Arm 1 (positive, transcription): p53 protein activates Mdm2 mRNA transcription  (R4).
Arm 2 (negative, protein):       Mdm2 binds p53 -> complex -> p53 degraded        (R8 + R10).
Nutlin-3: occupies the p53 pocket of MDM2 -> blocks R8 (k_bin -> k_bin*(1-nutlin)) -> cuts arm 2.
"""

import numpy as np

species = ('p53_mRNA', 'p53', 'Mdm2_mRNA', 'Mdm2', 'Mdm2_p53')
IDX = {s: i for i, s in enumerate(species)}

# ---------------------------------------------------------------------------
# Parameters (units: /min, molecule counts). Structure & ratios from Proctor & Gray,
# tuned toward clear sustained oscillations.
# ---------------------------------------------------------------------------
DEFAULT_PARAMS = dict(
    k_txn_p53     = 1.0,    # R0  p53 transcription (CONSTANT - p53 regulated post-translationally)
    k_deg_p53mRNA = 0.10,   # R1  p53 mRNA degradation
    k_tsl_p53     = 2.0,    # R2  p53 translation
    k_deg_p53     = 0.02,   # R3  basal p53 degradation (slow)
    k_txn_mdm2    = 0.030,  # R4  p53-activated Mdm2 transcription (per p53) -- arm 1
    k_txn_mdm2_0  = 0.01,   # R4b basal (leaky) Mdm2 transcription
    k_deg_mdm2mRNA= 0.08,   # R5  Mdm2 mRNA degradation
    k_tsl_mdm2    = 1.0,    # R6  Mdm2 translation
    k_deg_mdm2    = 0.08,   # R7  Mdm2 protein degradation
    k_bin         = 0.006,  # R8  p53 + Mdm2 -> complex (binding)   <-- NUTLIN blocks this
    k_rel         = 0.01,   # R9  complex -> p53 + Mdm2 (dissociation)
    k_degp53      = 0.50,   # R10 p53 degraded within the complex (Mdm2 recycled) -- arm 2
)

# Initial molecule counts (basal, unstressed state).
initial_state = np.array([10, 20, 5, 5, 0])

# ---------------------------------------------------------------------------
# Stoichiometry matrix: transitions[:, j] = state change when reaction j fires.
# Species order: (p53_mRNA, p53, Mdm2_mRNA, Mdm2, Mdm2_p53)
# ---------------------------------------------------------------------------
transitions = np.array([
    ( 1, 0, 0, 0, 0),  # R0  0 -> p53_mRNA                (p53 transcription)
    (-1, 0, 0, 0, 0),  # R1  p53_mRNA -> 0                (mRNA degradation)
    ( 0, 1, 0, 0, 0),  # R2  p53_mRNA -> p53_mRNA + p53   (translation)
    ( 0,-1, 0, 0, 0),  # R3  p53 -> 0                     (basal degradation)
    ( 0, 0, 1, 0, 0),  # R4  p53 -> p53 + Mdm2_mRNA       (p53-activated Mdm2 txn)  [arm 1]
    ( 0, 0, 1, 0, 0),  # R4b 0 -> Mdm2_mRNA               (basal Mdm2 transcription)
    ( 0, 0,-1, 0, 0),  # R5  Mdm2_mRNA -> 0               (mRNA degradation)
    ( 0, 0, 0, 1, 0),  # R6  Mdm2_mRNA -> Mdm2_mRNA+Mdm2  (translation)
    ( 0, 0, 0,-1, 0),  # R7  Mdm2 -> 0                    (protein degradation)
    ( 0,-1, 0,-1, 1),  # R8  p53 + Mdm2 -> complex        (binding)   [Nutlin blocks]
    ( 0, 1, 0, 1,-1),  # R9  complex -> p53 + Mdm2        (dissociation)
    ( 0, 0, 0, 1,-1),  # R10 complex -> Mdm2             (p53 degraded, Mdm2 recycled) [arm 2]
]).T

NUM_REACTIONS = transitions.shape[1]


def make_propensities(params=None, nutlin=0.0):
    """
    Return a propensities(x) function in the Patterns-2021 style.
    nutlin in [0,1]: Nutlin-3 efficacy; multiplies k_bin by (1-nutlin) (blocks p53-Mdm2 binding).
    """
    p = dict(DEFAULT_PARAMS)
    if params:
        p.update(params)
    kbin_eff = p['k_bin'] * (1.0 - nutlin)

    def propensities(x):
        return np.array([
            p['k_txn_p53'],                       # R0
            p['k_deg_p53mRNA'] * x[0],            # R1
            p['k_tsl_p53']     * x[0],            # R2
            p['k_deg_p53']     * x[1],            # R3
            p['k_txn_mdm2']    * x[1],            # R4  (p53-dependent)
            p['k_txn_mdm2_0'],                    # R4b
            p['k_deg_mdm2mRNA']* x[2],            # R5
            p['k_tsl_mdm2']    * x[2],            # R6
            p['k_deg_mdm2']    * x[3],            # R7
            kbin_eff           * x[1] * x[3],     # R8  (mass action p53*Mdm2)
            p['k_rel']         * x[4],            # R9
            p['k_degp53']      * x[4],            # R10
        ])
    return propensities


def ode_rhs(t, y, params=None, nutlin=0.0):
    """
    Mean-field (deterministic ODE) approximation of the same network -- for quickly scanning the
    oscillatory parameter region before running SSA. y = [p53_mRNA, p53, Mdm2_mRNA, Mdm2, complex].
    """
    p = dict(DEFAULT_PARAMS)
    if params:
        p.update(params)
    kbin_eff = p['k_bin'] * (1.0 - nutlin)
    mp, P, md, D, C = y
    bind = kbin_eff * P * D
    dmp = p['k_txn_p53'] - p['k_deg_p53mRNA'] * mp
    dP  = p['k_tsl_p53'] * mp - p['k_deg_p53'] * P - bind + p['k_rel'] * C
    dmd = p['k_txn_mdm2'] * P + p['k_txn_mdm2_0'] - p['k_deg_mdm2mRNA'] * md
    dD  = p['k_tsl_mdm2'] * md - p['k_deg_mdm2'] * D - bind + p['k_rel'] * C + p['k_degp53'] * C
    dC  = bind - p['k_rel'] * C - p['k_degp53'] * C
    return [dmp, dP, dmd, dD, dC]
