# Theory: the p53–MDM2 feedback loop, Nutlin-3, and the moment framework

## 1. The p53–MDM2 negative feedback loop

p53 ("guardian of the genome") is a transcription-factor tumour suppressor. Its abundance is
controlled by a two-arm loop with MDM2:

- **Arm 1 — transcriptional (positive):** p53 protein binds the *MDM2* promoter and **increases
  MDM2 transcription** ⇒ more MDM2 mRNA ⇒ more MDM2 protein.
- **Arm 2 — post-translational (negative):** MDM2 is an E3 ubiquitin ligase; it binds p53,
  ubiquitinates it, and targets it for proteasomal **degradation** ⇒ less p53.

Chained: p53 ↑ → MDM2 ↑ → p53 ↓ → MDM2 ↓ → p53 ↑ … = a **negative feedback loop**.

**Key facts used in the model**
- p53 is regulated mostly *post-translationally*: transcription of *p53* is ~constant; what sets
  p53 level is its **degradation rate**, controlled by MDM2. So in the model *p53 transcription is
  constant* while *MDM2 transcription depends on p53 protein*.
- The **MDM2 mRNA** step introduces a **time delay** (transcription must precede translation).
  This delay is sufficient for the loop to produce **sustained oscillations** after stress
  (Proctor & Gray 2008). Deterministically the loop is a *damped* oscillator; **intrinsic noise**
  in a stochastic (Markov-jump) description sustains the oscillations — matching single-cell p53
  "pulses" seen experimentally (Lahav, Geva-Zatorsky).

## 2. Nutlin-3 opens the loop

**Nutlin-3** (and the clinical analogue **Idasanutlin / RG7388**) is a small molecule that inserts
into the **p53-binding pocket of MDM2**. Consequences:
1. MDM2 can no longer bind p53 ⇒ arm 2 (degradation) is **cut**.
2. p53 **accumulates**.
3. Arm 1 is intact, so high p53 keeps driving MDM2 transcription ⇒ **MDM2 mRNA rises** (bounded by
   promoter saturation), but the MDM2 protein made **cannot restrain p53** ⇒ the loop is **open**,
   oscillations vanish.

Because Nutlin acts through p53, it activates p53 targets **only in TP53 wild-type cells**; in
TP53-mutant cells nothing happens (used as a biological negative control in step 7).

In the model, Nutlin efficacy `nutlin ∈ [0,1]` multiplies the MDM2→p53 degradation rate by
`(1 − nutlin)`.

## 3. Reference stochastic model (Proctor & Gray 2008)

A molecular-level Markov-jump model simulated with the Gillespie algorithm that reproduces
p53–MDM2 oscillations, using the MDM2 mRNA intermediate as the delay. Our `src/model.py` follows
this structure but uses the classic oscillator propensities (Michaelis–Menten degradation of p53
by MDM2, Hill activation of MDM2 transcription by p53) tuned to oscillate; `src/model_elementary.py`
keeps a fully elementary mass-action variant.

## 4. Why moments (the Patterns 2021 idea)

Only **mRNA** is observed in scRNA-seq. The regulatory signature is spread across the **statistical
moments** of the mRNA counts — mean (order 1), covariance (order 2), skewness (order 3) — not the
mean alone. The raw moment of a snapshot at time `t_d` is estimated as

    E[ A^{l1} B^{l2} ](t_d) ≈ (1/N) Σ_n a_{d,n}^{l1} b_{d,n}^{l2}

with A = p53 mRNA, B = MDM2 mRNA. Recovering *directed* regulation from these moments is what the
paper's non-linear **moment-based inference (MBI)** does (step 8).

## References
- Raharinirina et al. (2021) *Patterns* 2, 100332.
- Proctor & Gray (2008) *BMC Systems Biology* 2:75.
- McFarland et al. (2020) *Nat Commun* 11, 4296 (MIX-seq).
- Zappia et al. (2017) *Genome Biology* 18:174 (Splatter).
- Vassilev et al. (2004) *Science* 303:844 (Nutlin).
