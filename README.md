# Simulating the p53–MDM2 feedback loop (+ Nutlin-3) with the *Patterns 2021* moment framework

This project reuses the methodology of **Raharinirina et al., _Patterns_ 2021**
("Inferring gene regulatory networks from single-cell RNA-seq temporal snapshot data requires
higher-order moments") to **simulate the p53–MDM2 gene regulatory loop**, model the drug
**Nutlin-3 / Idasanutlin** (which opens the loop), and test a range of network-inference ideas
on both **synthetic** and **real** single-cell RNA-seq data.

> **Full write-up:** open [`report/report.html`](report/report.html) — a self-contained report
> (theory, methods with equations, and all results explained simply).

---

## What is inside

| Layer | Files |
|---|---|
| **Core library** | `src/ssa.py` (original Patterns2021 Gillespie engine), `src/model.py` (p53–MDM2 model), `src/model_elementary.py` (alternative mass-action model), `src/fast_ssa.py` (numba-accelerated ensemble simulator), `src/moments.py`, `src/association.py`, `src/directional.py` |
| **Vendored MBI** | `src/mbi/` — the original non-linear moment-based inference code (see `src/mbi/ATTRIBUTION.md`) |
| **Analysis steps** | `analysis/step1..step9` |
| **Outputs** | `figures/*.png`, `data/*.npz` |
| **Docs / report** | `docs/theory.md`, `report/report.html` |

## The biology in one paragraph

p53 is a tumour-suppressor transcription factor. It **activates transcription of MDM2**
(positive/transcriptional arm). MDM2 protein in turn **ubiquitinates p53 and sends it for
degradation** (negative/post-translational arm). Together these form a **delayed negative
feedback loop** (the delay comes from the MDM2 mRNA→protein step), which produces **pulses /
oscillations** of p53 after stress. **Nutlin-3** (and its clinical analogue **Idasanutlin**)
binds the p53-pocket of MDM2, blocks MDM2→p53 degradation, so p53 **accumulates**, MDM2 mRNA
**rises** (p53 still drives its transcription) but can no longer restrain p53 → the loop is
**open**. See `docs/theory.md` for details and equations.

## The steps

1. **step1** – single-cell dynamics (ODE damped vs SSA noise-sustained oscillations; closed vs Nutlin).
1b. **step1b** – validate `fast_ssa` ≡ `ssa.py` (identical propensities/stoichiometry).
2. **step2** – Nutlin dose response (p53↑, MDM2 mRNA↑, oscillation↓).
3. **step3** – higher-order moments (mean/covariance/skewness) of the (p53, MDM2) mRNA pair.
4. **step4** – Figure-3-style plots (snapshot scatter + moment time courses).
5. **step5** – add simple dropout noise, compare Pearson / MI / HSIC.
6. **step6** – **published Splatter noise** + Pearson/Spearman/MI/HSIC/**dCor** + permutation tests.
7. **step7 (a–d)** – **real scRNA-seq** (MIX-seq): load & QC (7a); DMSO vs Idasanutlin with negative
   controls + a WT/mutant 2×2 control (7b); 6 h vs 24 h timepoints (7c); the direct (TP53, MDM2) pair
   shown honestly to be weak on mRNA (7d).
8. **step8** – the non-linear **MBI** of Raharinirina et al. (2021) recovers the directed edge p53→MDM2 from mRNA moments.
9. **step9** – make Cor/MI/HSIC/dCor **directional** (lag, Granger, transfer entropy) and separate
   **regulation vs correlation** (partial correlation / conditioning). **step9b** saves a full
   six-species p53-MDM2-CDKN1A simulation under Nutlin for Figure 13D and later reuse.

## Install & run

```bash
pip install -r requirements.txt

# each step is self-contained; run from the repo root
python analysis/step1_single_cell_dynamics.py
python analysis/step3_scRNAseq_moments.py
python analysis/step8_mbi_inference.py
# ... etc. Figures land in figures/, data in data/.

# rebuild the HTML report after (re)generating figures
python report/build_report.py
```

Each step caches its (expensive) simulation to `data/cache_stepN.npz` and separates
`compute()` from `plot()`, so re-running only re-plots (seconds). Force a fresh simulation with
`RECOMPUTE=1 python analysis/stepN.py`. The engine is seeded per cell, so a given seed reproduces
the ensemble exactly.

### Real-data step (7)
`analysis/step7a_load_realdata.py` expects the **MIX-seq** dataset (McFarland et al.,
_Nat Commun_ 2020) in a folder `10298696/` next to this repo, containing the experiment
subfolders `DMSO_6hr_expt1/` and `Idasanutlin_6hr_expt1/` (10x `matrix.mtx` + `genes.tsv` +
`barcodes.tsv` + `classifications.csv`). A compact extract is cached in
`data/data_step7_panel.npz`, which `step7b` reads directly.

## Key results (short)

- Intrinsic noise turns the deterministically-damped loop into **sustained p53 pulses**; Nutlin
  makes p53 accumulate and abolishes oscillations.
- The regulatory signal lives in **higher-order moments** (covariance, skewness), not the mean.
- Under realistic (Splatter) noise, **library-size** creates spurious gene–gene correlations that
  fool all the powerful measures (Pearson/Spearman/HSIC/dCor); only the low-power kNN mutual-information
  estimator stays near zero. The real fix is **library-size normalization**, not the choice of statistic.
- On **real** MIX-seq data, MDM2–CDKN1A co-expression appears **only in TP53-WT cells under
  Idasanutlin** (dCor 0.12→0.38), with a clean mutant-line negative control. (The direct TP53–MDM2 pair
  is weak because TP53 mRNA is a poor proxy for p53 activity.)
- **Non-linear MBI recovers the directed edge p53→MDM2** from mRNA moments alone; this direction is
  **robust across random seeds**, though the exact edge weight is not (report direction, not magnitude).
- Symmetric measures become **directional** with time (lag/Granger/transfer entropy) and separate
  **regulation from correlation** with **conditioning** (partial correlation) — but these need
  per-cell time series, which real snapshots lack (hence moment dynamics / RNA velocity).

## References
- Raharinirina et al. (2021) *Patterns* 2, 100332. Code: https://github.com/vikramsunkara/ScRNAseqMoments
- Proctor & Gray (2008) *BMC Systems Biology* 2:75 (stochastic p53–Mdm2 oscillations).
- McFarland et al. (2020) *Nat Commun* 11, 4296 (MIX-seq).
- Zappia et al. (2017) *Genome Biology* 18:174 (Splatter).
- Vassilev et al. (2004) *Science* 303:844 (Nutlin).

See `LICENSE` (MIT for this project; `src/mbi/` keeps its upstream authors' rights).
