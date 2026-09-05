# Attribution: `src/mbi/` (moment-based inference)

The Python files in this `mbi/` package are **vendored (copied) from the original repository**
accompanying the paper:

> N. A. Raharinirina, F. Peppert, M. von Kleist, C. Schütte, V. Sunkara (2021).
> *Inferring gene regulatory networks from single-cell RNA-seq temporal snapshot data
> requires higher-order moments.* **Patterns** 2, 100332.
> https://doi.org/10.1016/j.patter.2021.100332

Original code: **https://github.com/vikramsunkara/ScRNAseqMoments**
(directory `MBI_Methods/Methods`).

Files included here: `two_Dim_mRNA.py`, `Symbolic_Moment_Generator.py`, `util.py`,
`Der_Based_Spline.py`, `NLLS.py`, `NLLS_Push_Forward.py`, `Regulatory_Rules.py`.

These files implement the **non-linear moment-based inference (MBI)** used in `step8`.
They are the intellectual property of the original authors and are redistributed here only to
make the analysis reproducible. Please refer to the upstream repository for its license and
cite the paper above if you use this method. The `linear MBI` (SINDy/cvxopt) path from the
original repo is **not** vendored here (step 8 uses a `scipy.optimize.nnls` initializer instead,
which needs no cvxopt).
