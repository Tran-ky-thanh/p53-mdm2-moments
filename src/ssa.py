####
# Gillespie Stochastic Simulation Algorithm (SSA) engine.
#
# Vendored VERBATIM from the repository of the Patterns 2021 paper:
#   Raharinirina et al. (2021), Patterns 2, 100332.
#   github.com/vikramsunkara/ScRNAseqMoments  ->  GRN_Models/SSA.py
#   ("Code For Constructing Realisations of Kurtz Processes. Author The big V.")
#
# Kept identical to ensure we use exactly their framework.
####

import numpy as np


def Time_To_Next_Reaction(lam):
    """Sample the waiting time to the next reaction ~ Exp(lam)."""
    r = np.random.rand()
    while r == 0:
        r = np.random.rand()
    return (1.0 / lam) * np.log(1.0 / r)


def Find_Reaction_Index(a):
    """Pick a reaction index from the cumulative propensities (roulette wheel)."""
    r = np.random.rand()
    while r == 0:
        r = np.random.rand()
    cum_prop = np.cumsum(a)
    j_chosen = np.sum(cum_prop < r * cum_prop[-1])
    return j_chosen


def SSA(Stochiometry, Propensities, X_0, t_0, t_final):
    """
    Gillespie algorithm. Given the stoichiometry, propensities and initial state, return a
    sample of the Kurtz process at t_final.

    Stochiometry : array (Num_species, Num_reaction).
    Propensities : function state -> propensity vector.
    X_0          : array (Num_species,).
    """
    t = t_0
    x = X_0.copy()

    while t <= t_final:
        a = Propensities(x)
        tau = Time_To_Next_Reaction(np.sum(a))
        # If we would jump past t_final, or total propensity is 0, stop and keep the state
        if (t + tau > t_final) or (np.sum(a) == 0):
            return x, t_final
        else:
            t = t + tau
            j = Find_Reaction_Index(a)
            x = x + Stochiometry[:, j]


def SSA_Fixed_Width_Trajectory(Stochiometry, Propensities, X_0, T_Obs_Points):
    """Sample a trajectory at fixed observation times T_Obs_Points -> snapshots."""
    X = [X_0]
    for i in range(len(T_Obs_Points) - 1):
        x, t = SSA(Stochiometry, Propensities, X[-1], T_Obs_Points[i], T_Obs_Points[i + 1])
        X.append(x)
    return np.array(X).T  # shape (Num_species, Num_time_points)
