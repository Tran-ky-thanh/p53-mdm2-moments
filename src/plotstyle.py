# -*- coding: utf-8 -*-
"""Shared matplotlib style: larger, readable fonts for all figures."""
import matplotlib as mpl


def apply():
    mpl.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 15,
        'axes.labelsize': 14,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'figure.titlesize': 17,
    })
