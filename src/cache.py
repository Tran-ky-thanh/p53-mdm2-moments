# -*- coding: utf-8 -*-
"""
Tiny caching helper so plotting does NOT re-run the (expensive) simulation every time.

Convention: each analysis step splits into
    compute() -> returns a flat dict of numpy arrays  (the expensive simulation)
    plot(D)   -> draws the figure from that dict       (cheap)
and calls  D = load_or_compute(cache_path, compute).

The first run simulates and writes data/cache_stepX.npz; later runs (e.g. when you only
tweak figure labels) just LOAD the cache and re-plot in seconds.

Force a recompute with  RECOMPUTE=1  in the environment, or `force=True`.
"""
import os
import pathlib
import numpy as np


def load_or_compute(path, compute_fn, force=None):
    path = pathlib.Path(path)
    if force is None:
        force = os.environ.get("RECOMPUTE", "0") == "1"
    if path.exists() and not force:
        with np.load(path, allow_pickle=True) as d:
            data = {k: d[k] for k in d.files}
        print(f"[cache] loaded {path.name} ({len(data)} arrays) - skipping simulation")
        return data
    print(f"[cache] computing {path.name} ...")
    data = compute_fn()
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **data)
    print(f"[cache] saved {path.name}")
    return data
