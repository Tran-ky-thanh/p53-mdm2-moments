# -*- coding: utf-8 -*-
"""Six-species p53-MDM2-CDKN1A extension used by analysis step 9.

MDM2 and CDKN1A are both transcriptional targets of p53 protein. MDM2
retains its negative-feedback action on p53; CDKN1A is a downstream target
and does not feed back in this minimal model.
"""
import numpy as np
import model as BASE

species = ('p53_mRNA', 'p53', 'Mdm2_mRNA', 'Mdm2', 'CDKN1A_mRNA', 'CDKN1A')
IDX = {s: i for i, s in enumerate(species)}

DEFAULT_PARAMS = dict(BASE.DEFAULT_PARAMS)
DEFAULT_PARAMS.update(
    k_txn_cdkn1a=5.0, Kd_cdkn1a=70.0, n_cdkn1a=3.0,
    k_txn_cdkn1a_0=0.05, k_deg_cdkn1a_mRNA=0.08,
    k_tsl_cdkn1a=1.0, k_deg_cdkn1a=0.15,
)

initial_state = np.array([15, 30, 20, 50, 15, 40])

# R0-R8 match model.py; R9-R12 add CDKN1A transcription, mRNA decay,
# translation and protein decay.
transitions = np.array([
    ( 1, 0, 0, 0, 0, 0), (-1, 0, 0, 0, 0, 0),
    ( 0, 1, 0, 0, 0, 0), ( 0,-1, 0, 0, 0, 0),
    ( 0,-1, 0, 0, 0, 0), ( 0, 0, 1, 0, 0, 0),
    ( 0, 0,-1, 0, 0, 0), ( 0, 0, 0, 1, 0, 0),
    ( 0, 0, 0,-1, 0, 0), ( 0, 0, 0, 0, 1, 0),
    ( 0, 0, 0, 0,-1, 0), ( 0, 0, 0, 0, 0, 1),
    ( 0, 0, 0, 0, 0,-1),
]).T


def make_propensities(params=None, nutlin=0.0):
    """Return exact SSA hazards for the six-species extension."""
    p = dict(DEFAULT_PARAMS)
    if params:
        p.update(params)
    kdeg_eff = p['k_degp53'] * (1.0 - nutlin)

    def propensities(x):
        mp, P, md, D, mc, C = x
        hp = P**p['n_hill'] / (p['Kd_mdm2']**p['n_hill'] + P**p['n_hill'])
        hc = P**p['n_cdkn1a'] / (p['Kd_cdkn1a']**p['n_cdkn1a'] + P**p['n_cdkn1a'])
        return np.array([
            p['k_txn_p53'], p['k_deg_p53mRNA'] * mp, p['k_tsl_p53'] * mp,
            p['k_deg_p53'] * P, kdeg_eff * D * P / (p['Km_p53'] + P),
            p['k_txn_mdm2'] * hp + p['k_txn_mdm2_0'],
            p['k_deg_mdm2mRNA'] * md, p['k_tsl_mdm2'] * md, p['k_deg_mdm2'] * D,
            p['k_txn_cdkn1a'] * hc + p['k_txn_cdkn1a_0'],
            p['k_deg_cdkn1a_mRNA'] * mc, p['k_tsl_cdkn1a'] * mc,
            p['k_deg_cdkn1a'] * C,
        ])
    return propensities
