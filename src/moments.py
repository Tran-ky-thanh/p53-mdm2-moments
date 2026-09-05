# -*- coding: utf-8 -*-
"""
Utilities to compute MOMENTS of mRNA counts over time -- in the spirit of Raharinirina et al. (2021).

Following Equation 1 of the paper, the raw moment of order l = l1 + l2 at snapshot d is
    E[ A^l1 B^l2 ](t_d) ~ (1/N) * sum_n a_{d,n}^l1 * b_{d,n}^l2
with A, B the counts of the two mRNAs (here A = p53 mRNA, B = MDM2 mRNA) and N the cell number.

Provides: order-1 mean, order-2 variance/covariance, order-3 standardized skewness
(mean + covariance + skewness together carry the regulatory information, per Raharinirina et al. 2021).
"""
import numpy as np


def raw_moment(A, B, l1, l2):
    """E[A^l1 B^l2] at each snapshot. A, B have shape (N_cells, N_time)."""
    return np.mean((A ** l1) * (B ** l2), axis=0)


def moment_time_courses(A, B):
    """
    A, B: arrays (N_cells, N_time) of mRNA counts (p53, MDM2) per cell and time.
    Return a dict of statistic time-courses.
    """
    N = A.shape[0]
    mA = A.mean(0); mB = B.mean(0)                         # order-1 mean
    varA = A.var(0); varB = B.var(0)                       # order-2 variance
    cov = ((A - mA) * (B - mB)).mean(0)                    # order-2 covariance
    sdA = np.sqrt(varA); sdB = np.sqrt(varB)
    with np.errstate(divide='ignore', invalid='ignore'):
        skA = (((A - mA) ** 3).mean(0)) / (sdA ** 3)      # order-3 standardized skewness
        skB = (((B - mB) ** 3).mean(0)) / (sdB ** 3)
        corr = cov / (sdA * sdB)
    return dict(mean_p53=mA, mean_mdm2=mB, var_p53=varA, var_mdm2=varB,
                cov=cov, corr=corr, skew_p53=skA, skew_mdm2=skB, N=N)
