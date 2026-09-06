# -*- coding: utf-8 -*-
"""
Build a self-contained report/report.html from the figures in figures/.
Figures are embedded as base64 so the HTML works offline and on GitHub Pages.

Run:  python report/build_report.py
"""
import base64, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
OUT = ROOT / "report" / "report.html"


def img(name):
    p = FIG / name
    if not p.exists():
        return f'<div class="missing">[missing figure: {name} - run the corresponding step]</div>'
    b64 = base64.b64encode(p.read_bytes()).decode()
    return f'<img src="data:image/png;base64,{b64}" alt="{name}">'


SCHEMATIC = '''
<svg viewBox="0 0 720 260" width="100%" style="max-width:720px">
  <defs>
    <marker id="arr" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto">
      <path d="M0,0 L8,3 L0,6 Z" fill="#1a6"/></marker>
    <marker id="bar" markerWidth="12" markerHeight="12" refX="1" refY="6" orient="auto">
      <path d="M1,1 L1,11" stroke="#c33" stroke-width="3"/></marker>
  </defs>
  <rect x="70" y="100" width="150" height="60" rx="10" fill="#dbeafe" stroke="#2563eb"/>
  <text x="145" y="135" text-anchor="middle" font-size="20" font-weight="bold" fill="#1e3a8a">p53</text>
  <rect x="500" y="100" width="150" height="60" rx="10" fill="#fde2e4" stroke="#be123c"/>
  <text x="575" y="135" text-anchor="middle" font-size="20" font-weight="bold" fill="#9f1239">MDM2</text>
  <path d="M220,115 C330,60 430,60 500,115" fill="none" stroke="#1a6" stroke-width="3" marker-end="url(#arr)"/>
  <text x="360" y="55" text-anchor="middle" font-size="14" fill="#1a6">Arm 1: transcription (+)  p53 &#8594; MDM2 mRNA</text>
  <path d="M500,150 C430,205 330,205 224,150" fill="none" stroke="#c33" stroke-width="3" marker-end="url(#bar)"/>
  <text x="360" y="225" text-anchor="middle" font-size="14" fill="#c33">Arm 2: degradation (&#8722;)  MDM2 &#8867; p53</text>
  <circle cx="360" cy="198" r="16" fill="#fff7cc" stroke="#eab308" stroke-width="2"/>
  <text x="360" y="203" text-anchor="middle" font-size="11" font-weight="bold" fill="#a16207">Nutlin</text>
  <line x1="345" y1="183" x2="375" y2="213" stroke="#eab308" stroke-width="3"/>
</svg>'''

CSS = '''
:root{--fg:#1f2937;--muted:#6b7280;--bg:#ffffff;--card:#f8fafc;--accent:#2563eb;--line:#e5e7eb}
*{box-sizing:border-box}
body{margin:0;font:16px/1.7 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:var(--fg);background:#eef2f7}
.wrap{max-width:1000px;margin:0 auto;background:var(--bg);padding:0 0 60px}
header{background:linear-gradient(135deg,#1e3a8a,#2563eb);color:#fff;padding:40px 48px}
header h1{margin:0 0 8px;font-size:27px;line-height:1.3}
header p{margin:0;opacity:.92}
main{padding:0 48px}
h2{margin-top:46px;padding-top:14px;border-top:2px solid var(--line);font-size:23px;color:#0f172a}
h3{margin-top:34px;font-size:19px;color:#1e3a8a;line-height:1.35}
p,li{color:var(--fg)}
.lead{font-size:17px;color:#374151}
.toc{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px 22px;margin:24px 0;font-size:14.5px}
.toc a{color:var(--accent);text-decoration:none}.toc a:hover{text-decoration:underline}
.formula{background:#0f172a;color:#e2e8f0;border-radius:8px;padding:12px 16px;margin:14px 0;
  font-family:"SF Mono",Consolas,monospace;font-size:14.5px;overflow-x:auto;white-space:pre}
.note{background:#fff7ed;border-left:4px solid #f59e0b;padding:10px 16px;border-radius:6px;margin:14px 0}
.key{background:#ecfdf5;border-left:4px solid #10b981;padding:10px 16px;border-radius:6px;margin:14px 0}
.samples{background:#eef2ff;border:1px solid #c7d2fe;border-radius:6px;padding:7px 13px;margin:12px 0;
  font-size:13.5px;color:#3730a3}
figure{margin:22px 0 6px}
figure img{width:100%;border:1px solid var(--line);border-radius:8px;display:block;background:#fff}
figcaption{font-size:14.5px;color:#374151;margin-top:10px;padding:0 4px}
figcaption .lab{font-weight:700;color:#0f172a}
table{border-collapse:collapse;width:100%;margin:14px 0;font-size:14.5px}
th,td{border:1px solid var(--line);padding:8px 10px;text-align:left}
th{background:var(--card)}
.svgbox{text-align:center;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:18px;margin:18px 0}
code{background:#eef2f7;padding:1px 5px;border-radius:4px;font-size:14px}
.missing{color:#b91c1c;background:#fee2e2;padding:10px;border-radius:6px}
.refs{font-size:14.5px}
footer{color:var(--muted);font-size:13.5px;padding:30px 48px;border-top:1px solid var(--line)}
'''

# ---- Results: each entry = academic subsection with caption + interpretation ----
RESULTS = [
 dict(sid="r1", num="4.1", fig="fig_step1_single_cell.png", fign=1,
   title="Stochastic vs deterministic dynamics: intrinsic noise sustains p53 pulses, and Nutlin-3 opens the loop",
   samples="Deterministic ODE integration plus <b>4 stochastic (Gillespie) cells per condition</b>; horizon 1500 min, sampled every 2 min. Few single cells are shown on purpose so individual traces remain legible.",
   caption="Deterministic (top) and stochastic (bottom) trajectories of the p53-MDM2 circuit in the closed loop (left) versus under full Nutlin-3 inhibition (right). Blue: p53 protein; orange: MDM2 mRNA; green: MDM2 protein.",
   interp=[
    "In the closed loop the mean-field ODE exhibits a <i>damped</i> oscillation that relaxes to a stable focus, whereas the stochastic simulation of the identical reaction network produces <i>sustained</i>, irregular oscillations. This dichotomy&mdash;deterministically damped but stochastically sustained&mdash;is the mechanism by which discrete, noisy molecular events continually re-excite the loop, and it reproduces the discrete p53 'pulses' reported in single living cells (Lahav et al., 2004; Geva-Zatorsky et al., 2006). It is the same qualitative behaviour identified in the mechanistic stochastic model of Proctor &amp; Gray (2008).",
    "Full Nutlin-3 inhibition (right) removes the MDM2-mediated degradation of p53 (arm 2). p53 no longer has a fast decay channel and rises to a high plateau in every cell, while the oscillation disappears: the negative feedback is broken and the system is driven to an open-loop steady state. This is the single-cell signature of pharmacological MDM2 inhibition (Vassilev et al., 2004)."]),
 dict(sid="r2", num="4.2", fig="fig_step2_dose_response.png", fign=2,
   title="Dose-dependent p53 accumulation, MDM2-mRNA induction, and suppression of oscillation",
   samples="<b>1500 cells per dose</b>, six Nutlin-3 efficacies; horizon 1600 min; statistics in the late window t &ge; 600 min. "
           "Here 'efficacy' is a mechanistic parameter (fraction of p53-MDM2 binding blocked, 0-1), <i>not</i> a drug concentration; "
           "the experimental data (Section 4.7) used a single near-saturating dose of 2.5 uM (no dose series was available).",
   caption="Population statistics in the late window as a function of Nutlin-3 efficacy (0 = vehicle, 1 = complete block of p53-MDM2 binding): mean p53 protein (left), mean MDM2 mRNA (centre), and a single-cell oscillation index defined as the mean temporal coefficient of variation of p53 (right).",
   interp=[
    "Increasing inhibition raises mean p53 monotonically and steeply near full block (left), consistent with progressive loss of the degradation channel. MDM2 mRNA also rises (centre) but saturates, because transcriptional activation of MDM2 by p53 (arm 1) is bounded by promoter saturation&mdash;p53 keeps driving transcription even though the resulting MDM2 protein can no longer act on p53.",
    "The oscillation index (right) is essentially flat at low-to-moderate doses and then collapses near complete inhibition, marking the transition from a functioning (pulsatile) loop to an open loop. Together the three curves quantify the qualitative picture of Fig. 1: arm 1 intact, arm 2 cut."]),
 dict(sid="r3", num="4.3", fig="fig_step3_moments.png", fign=3,
   title="Regulatory information resides in the higher-order moments of mRNA counts",
   samples="<b>3000 cells per condition</b>; snapshots every 15 min up to 1500 min. Only mRNA counts are used, mimicking scRNA-seq.",
   caption="Time courses of the first three statistical moments of the (p53 mRNA, MDM2 mRNA) pair across the cell population, in the closed loop versus under Nutlin-3: means (top), variance and covariance (middle), and skewness (bottom).",
   interp=[
    "The mean of p53 mRNA is essentially <i>identical</i> in the two conditions (the two lines overlap near 15 counts), because Nutlin acts on p53 protein stability, not on p53 transcription. An analysis restricted to the mean of this transcript would therefore detect no effect.",
    "The regulatory change is instead written into the higher-order moments: the covariance of the two mRNAs drops from about +10 (closed loop) to approximately 0 (Nutlin), and the skewness of MDM2 mRNA falls markedly. This is a direct, quantitative illustration of the central thesis of Raharinirina et al. (2021): the signature of regulation is distributed across mean, covariance and skewness rather than being captured by the mean alone."]),
 dict(sid="r4", num="4.4", fig="fig_step4_figure3_style.png", fign=4,
   title="Joint distribution and moment dynamics of the (p53, MDM2) mRNA pair",
   samples="<b>4000 cells per condition</b>; the scatter panels show a snapshot at t = 300 min.",
   caption="Single-cell joint distribution of p53 and MDM2 mRNA with density contours and the population mean trajectory (orange), for the closed loop (A) and Nutlin-3 (B); mean counts (C), covariance (D) and skewness (E, F) over time. This presentation parallels the moment-based summaries of Raharinirina et al. (2021).",
   interp=[
    "In the closed loop the joint cloud is tilted along a positive diagonal (covariance +10.5, correlation +0.43): cells rich in p53 mRNA also tend to be rich in MDM2 mRNA, the fingerprint of an intact transcriptional coupling. Under Nutlin the cloud becomes round and shifts upward (covariance approximately 0, correlation approximately 0), with MDM2 mRNA elevated&mdash;dependence is lost even though the marginal level of MDM2 mRNA is higher.",
    "Panel C again shows that the two solid p53-mRNA curves coincide, underscoring that the discriminating information lies in the second- and third-order statistics (D-F), not in the first-order mean of p53 mRNA."]),
 dict(sid="r5", num="4.5", fig="fig_step5_noise_association.png", fign=5,
   title="Detecting statistical dependence under dropout: correlation, mutual information and HSIC",
   samples="<b>4000 cells per condition</b>; snapshot at t = 600 min; five capture efficiencies &times; six independent noise draws each.",
   caption="Dependence between p53 and MDM2 mRNA estimated by Pearson correlation, mutual information (kNN estimator) and the Hilbert-Schmidt Independence Criterion (HSIC), as a function of a simple binomial capture (dropout) probability, in the closed loop versus Nutlin-3.",
   interp=[
    "All three measures report clear dependence in the closed loop and values near zero under Nutlin-3, and all decay as dropout increases. In this dependence, which is approximately linear and positive, Pearson correlation gives the cleanest separation; mutual information and HSIC confirm the same conclusion.",
    "This step uses a deliberately idealized noise model (each mRNA molecule retained independently with probability &beta;). The following step replaces it with a published, more realistic model and reveals a confounder that this simple model does not expose."]),
 dict(sid="r6", num="4.6", fig="fig_step6_splatter_methods.png", fign=6,
   title="Realistic technical noise (Splatter) and the library-size confound",
   samples="<b>3000 cells per condition</b>; snapshot at t = 600 min; five dropout midpoints &times; five noise draws; kernel measures subsampled to 1000 cells.",
   caption="The same detection task under the Splatter technical-noise model (Zappia et al., 2017)&mdash;log-normal library size, negative-binomial overdispersion, and expression-dependent logistic dropout&mdash;evaluated with five dependence measures including distance correlation (dCor).",
   interp=[
    "Unlike the idealized model of Fig. 5, the Splatter model applies a shared per-cell library-size factor that multiplies every gene. This induces a <i>spurious</i> positive association between transcripts that are in fact independent, so the Nutlin-3 curves no longer sit at zero for the correlation-type measures. This is a well-known scRNA-seq confounder and a concrete example of the pre-processing pitfalls emphasized by Raharinirina et al. (2021).",
    "The permutation test in Fig. 7 quantifies which apparent detections are real and which are artefacts of library size."]),
 dict(sid="r7", num=None, fig="fig_step6_permutation_z.png", fign=7,
   title=None, samples=None,
   caption="Permutation z-scores (Section 3.8) for each measure at a moderate dropout level (N = 3000 cells; n = 2000 subsample; 200 permutations; averaged over 5 noise draws). Bars above the dashed line (z = 2, approximately p &lt; 0.05) count as detections.",
   interp=[
    "In the closed loop the true dependence is strongly significant for Pearson (z &approx; 6), Spearman (&approx; 5.5), HSIC (&approx; 12) and distance correlation (&approx; 12); mutual information is the exception (z &approx; 0.7), because the kNN MI estimator is a weak, high-variance detector on discrete count data (Section 4.5).",
    "Under Nutlin-3 the two transcripts are genuinely independent, so every bar there <i>should</i> be zero. Instead Pearson (z &approx; 3.5), Spearman (&approx; 2.9), HSIC (&approx; 3.6) and dCor (&approx; 5.3) are all significant <i>false positives</i>, produced by the shared per-cell library-size factor; only mutual information stays near zero (z &approx; 0.7). The correct reading is therefore not that any one statistic is intrinsically 'robust' &mdash; the powerful measures all detect the library-size artefact, and MI only escapes because of its low power. The real remedy is explicit <b>library-size normalization</b>, which we apply to the experimental data in Section 4.7."]),
 dict(sid="r8", num="4.7", title="Validation on experimental scRNA-seq: MIX-seq DMSO versus idasanutlin",
   blocks=[
    dict(sub="Why not the direct (p53, MDM2) pair, and why does the sign look reversed versus the simulation?",
      fig="fig_step7d_direct_pair.png", fign=8,
      samples="WT line LNCaP, 6 h (n=152 treated). dCor of TP53-MDM2 vs MDM2-CDKN1A.",
      caption="For the wild-type line: (A) mean expression of TP53, MDM2 and CDKN1A in DMSO vs idasanutlin; "
              "(B) scatter of the direct loop pair TP53 mRNA vs MDM2 mRNA; (C) dCor of TP53-MDM2 vs the "
              "p53-target pair MDM2-CDKN1A in each condition.",
      interp=[
       "Our ideal observable is the loop pair (p53, MDM2), and the simulation uses exactly that. On real scRNA-seq, however, <i>TP53 mRNA</i> is a poor readout of p53 <i>activity</i>: p53 is controlled post-translationally, so its transcript barely moves and does not track the protein that actually drives MDM2 transcription. Panel (A) confirms this directly&mdash;TP53 mRNA barely changes with the drug (it even dips slightly), while MDM2 and CDKN1A are strongly induced. Panels (B)-(C) show the direct TP53-MDM2 scatter has little structure (dCor 0.07 in DMSO, 0.10 under idasanutlin) and stays weak in both conditions, whereas the p53-<i>target</i> pair MDM2-CDKN1A jumps from 0.12 to 0.38. On mRNA-only data the loop's regulatory coherence is therefore best read from p53-target co-expression, not from the nominal (p53, MDM2) pair&mdash;an intrinsic limitation of measuring a post-translationally regulated hub at the mRNA level. We use MDM2-CDKN1A for the rest of this section.",
       "The <i>direction</i> of the drug effect on that pair also differs from the simulation, for a principled reason: cell-to-cell covariance requires both an active regulatory link <b>and</b> heterogeneity in the upstream driver. In the simulation the baseline (no drug) is an actively cycling loop, so covariance is high at baseline and Nutlin&mdash;by clamping p53&mdash;removes it. In the experiment the baseline (DMSO) is quiescent (low, uniform p53 activity), so target co-expression is low, and it is the drug that switches p53 on heterogeneously, creating the co-expression. Same principle, opposite baseline."]),
    dict(sub=None, fig="fig_step7_realdata.png", fign=9,
      samples="Experimental data from MIX-seq (McFarland et al., 2020): <b>1886 (DMSO) and 2634 (idasanutlin) QC-passing singlets</b> across a pool of 24 cell lines; the TP53 wild-type line LNCaP contributes n &approx; 100 / 152 cells; permutation tests use B = 200-300.",
      caption="Application to real scRNA-seq: DMSO (closed loop) versus idasanutlin, a clinical MDM2 inhibitor of the Nutlin family (open loop), after CP10k + log normalization. (A) induction of p53 target genes in a TP53-wild-type line but not a mutant line; (B, C) MDM2-CDKN1A joint distributions in the wild-type and mutant lines; (D) all five measures on the wild-type MDM2-CDKN1A pair; (E) the 2x2 wild-type/mutant biological control; (F) change in dCor upon treatment across gene pairs.",
      interp=[
       "Idasanutlin sharply induces canonical p53 targets (MDM2, CDKN1A, FDXR, ...) in the TP53-wild-type line LNCaP but not in the TP53-mutant line, as expected for a p53-dependent drug (McFarland et al., 2020). Concomitantly, MDM2 and CDKN1A&mdash;two direct p53 targets&mdash;become co-expressed under treatment: dCor rises from 0.12 (not significant) to 0.38 (permutation z &approx; 9). The 2x2 control (E) isolates the p53-dependent effect: only the wild-type + idasanutlin condition shows strong co-regulation, while the drug-treated mutant line does not.",
       "A caveat is made explicit in panel F: among all gene pairs, only MDM2-CDKN1A gains dependence upon treatment, whereas housekeeping pairs actually lose covariation&mdash;idasanutlin arrests the cell cycle in wild-type cells, removing a global source of gene-gene covariance. Housekeeping pairs are therefore <i>not</i> a perfect negative control in single-cell data; the cleanest control is the isogenic comparison against a non-responding (mutant) line."]),
    dict(sub="Does the effect generalize across cell lines?",
      fig="fig_step7e_across_lines.png", fign=10,
      samples="All <b>22 MIX-seq cell lines with &ge; 40 QC-passing cells in both DMSO and idasanutlin</b> (6 h). "
              "For each line: p53 response = mean induction of MDM2 + CDKN1A (idasanutlin minus DMSO, log-norm); "
              "dCor gain = distance correlation of (MDM2, CDKN1A) in idasanutlin minus in DMSO.",
      caption="One point per cell line: p53 response (x-axis) versus the gain in MDM2-CDKN1A co-expression "
              "under idasanutlin (y-axis), with a linear trend fitted across all 22 lines.",
      interp=[
       "There is a positive association (r = +0.44): lines whose p53 targets are induced more strongly by idasanutlin also tend to gain more MDM2-CDKN1A co-expression, consistent with the mechanism above. LNCaP, the wild-type line used throughout this section, sits at the extreme of both axes and is the cleanest example in the panel; CCFSTTG1, NCIH226 and DKMG (moderate p53 response) also show a modest positive gain.",
       "The relationship is noisy, however, and we report it as such rather than overstating it. Two lines with almost no p53 response (RCM1, BT549) show large <i>negative</i> dCor changes, most likely small-sample or line-specific noise rather than a p53-related effect, and they visibly pull down the correlation. With only 22 lines and per-line cell counts of a few hundred, this trend should be read as supportive but not as strong independent confirmation; it motivates, rather than replaces, the single-line analysis with its 2x2 biological control."]),
    dict(sub="Timepoint dependence: 6 h versus 24 h of idasanutlin",
      fig="fig_step7c_timepoints.png", fign=11,
      samples="Idasanutlin was applied at a single saturating concentration (<b>2.5 uM</b>) and assayed at "
              "<b>6 h and 24 h</b> (MIX-seq has no idasanutlin dose series). WT line LNCaP: 152 vs 102 treated cells.",
      caption="For the TP53 wild-type line, (A) induction of p53 target genes (idasanutlin minus DMSO) at "
              "6 h versus 24 h, and (B) MDM2-CDKN1A co-expression (dCor) at the two timepoints.",
      interp=[
       "The two readouts move in <i>opposite</i> directions with time. The <b>mean induction</b> of p53 targets is generally <b>larger at 24 h</b> (summed log-fold induction 6.6 at 6 h vs 9.4 at 24 h), reflecting continued accumulation of target mRNA. In contrast, the <b>cell-to-cell co-expression</b> of MDM2 and CDKN1A is <b>stronger at 6 h</b> (dCor 0.38) than at 24 h (0.26).",
       "This is expected: at 6 h the response is still in its heterogeneous transient, so p53 activity&mdash;and hence its targets&mdash;varies markedly from cell to cell, maximizing the covariation our measures detect; by 24 h the response has largely saturated and homogenized (and fewer cells are recovered, consistent with arrest/apoptosis), so the mean is higher but the covariation is lower. The rest of this section uses the 6 h timepoint, which gives the cleaner co-regulation signal. Note that 'efficacy' in our simulation (Section 4.2) is a mechanistic fraction-of-binding-blocked in [0,1], whereas 2.5 uM is a single, near-saturating <i>concentration</i>; the two are related monotonically but are not the same quantity."]),
   ]),
 dict(sid="r9", num="4.8", fig="fig_step8_mbi.png", fign=12,
   title="Recovering directed regulation from moments: non-linear moment-based inference (Raharinirina et al., 2021)",
   samples="<b>8000 cells per condition</b>; raw moments up to order 3-4; approximately 49 snapshots used for the fit after discarding the transient and subsampling. "
           "<b>Note on noise:</b> the input is <i>clean</i> SSA counts&mdash;they carry the <i>intrinsic</i> (molecular) noise that MBI actually exploits, but <i>no</i> technical scRNA-seq noise (dropout, library size) is applied here. This is a best-case test, as in Raharinirina et al. (2021).",
   caption="Non-linear moment-based inference (MBI) of Raharinirina et al. (2021), applied to the simulated moment time-courses of the (p53 mRNA, MDM2 mRNA) pair (original code, vendored in src/mbi). Top: model fit to the mean dynamics; middle: the inferred 2x2 regulatory network; bottom: inferred edge strengths.",
   interp=[
    "Unlike the symmetric dependence measures of the previous steps, MBI returns a <i>directed</i> network. In the closed loop it cleanly recovers the transcriptional edge p53 &#8594; MDM2 (A&#8594;B = +3.9) and essentially no MDM2 &#8594; p53 edge (B&#8594;A = 0). The absence of an MDM2 &#8594; p53 edge is expected, because the negative arm operates post-translationally on p53 protein and is invisible at the mRNA level. This demonstrates, on our system, the ability of MBI to infer regulatory direction from mRNA moments alone (Raharinirina et al., 2021).",
    "Under Nutlin-3 MBI still recovers a positive p53 &#8594; MDM2 edge (A&#8594;B = +2.0). This is biologically correct: Nutlin blocks the MDM2&ndash;p53 <i>protein</i> interaction, not the transcriptional arm, so p53 continues to drive MDM2 transcription (indeed MDM2 mRNA is elevated). However, the inference is now less clean&mdash;it also reports a spurious MDM2 &#8594; p53 edge (B&#8594;A = -1.1). The reason is identifiability: with the loop open, p53 is clamped at a high, near-constant level and MDM2 mRNA sits on a flat plateau, so the moment time-courses carry little <i>dynamic</i> information for the fit to exploit. MBI is therefore most reliable in the informative, oscillating closed-loop regime; a near-static system yields weaker, noisier network estimates.",
    "Finally, a scope note: MBI works on the moment time-courses and implicitly assumes the <i>observed</i> moments equal the <i>true</i> moments. It de-noises only by averaging over many cells (and by spline-smoothing the higher moments), which suppresses sampling noise&mdash;but it does not model technical scRNA-seq artefacts. Because dropout and library-size shifts distort exactly the higher-order moments MBI relies on (Section 4.6), such technical noise would bias the inference; the clean, intrinsic-noise-only setting here is deliberately a best case."]),
 dict(sid="r10", num="4.9", fig="fig_step9_directional.png", fign=13,
   title="Endowing dependence measures with direction and conditional specificity",
   samples="<b>4000 cells per condition</b> for the directional analysis (fluctuation window t &ge; 200 min). Panel D uses a separate <b>4000-cell, six-species p53-MDM2-CDKN1A Gillespie simulation</b> at 4-min snapshots from 0-80 min with Nutlin efficacy = 1; its complete output is saved in <code>data/cache_step9_three_gene_nutlin.npz</code>. Real-data partial correlation uses LNCaP (n &approx; 152).",
   caption="(A) lagged cross-correlation and (B) lagged distance correlation of the detrended p53/MDM2 mRNA "
           "fluctuations, shown for tau &ge; 0 as two directional curves; at positive lag the p53&rarr;MDM2 curve "
           "dominates and peaks near tau &asymp; +4 min while MDM2&rarr;p53 decays, i.e. p53 fluctuations precede "
           "MDM2. (C) Granger causality and transfer entropy; (D) marginal and partial correlations from an explicit "
           "p53-MDM2-CDKN1A mRNA/protein simulation under full Nutlin; (E) partial correlation of MDM2-CDKN1A "
           "given a p53-activity proxy on real data.",
   interp=[
    "Two ingredients convert the symmetric measures of Sections 4.5-4.6 into causal ones. <b>Direction requires time</b>: using per-cell trajectories, the lagged cross-correlation and lagged distance correlation peak at positive lag for p53&#8594;MDM2 (A leads B), and Granger causality and transfer entropy are strongly asymmetric (Granger A&#8594;B = 0.21 vs B&#8594;A = 0.00; transfer entropy 0.09 vs 0.001 nats), recovering the direction p53&#8594;MDM2. (Because these measures are symmetric, C<sub>B&#8594;A</sub>(&tau;)=C<sub>A&#8594;B</sub>(&minus;&tau;), it is enough to show <i>positive</i> lags only: over &tau; &ge; 0 the two directional curves are genuinely distinct, and the p53&#8594;MDM2 curve exceeding MDM2&#8594;p53 there means p53 fluctuations <i>precede</i> MDM2&mdash;the temporal-precedence signature of p53&#8594;MDM2. The ~4-min peak reflects the transcription/translation delay.) A further advantage of the <i>signed</i> cross-correlation (panel A) is the <b>negative lobe near &tau; &asymp; 16 min</b>: after p53 drives MDM2 up, the negative-feedback arm pulls the system back, giving an anti-correlation (overshoot) at roughly half the loop's oscillation timescale. This lobe is <i>consistent with</i> the closed negative-feedback loop&mdash;it is absent under Nutlin (open loop; grey dashed curve) and, being a change of <i>sign</i>, cannot appear in the non-negative distance correlation of panel B. <b>Panel D replaces the abstract common-driver cartoon with an explicit p53-MDM2-CDKN1A simulation.</b> It contains mRNA and protein for all three genes; full Nutlin cuts MDM2-mediated p53 degradation, while p53 protein continues to drive both target transcripts. All three correlations remain near zero during 0-80 min. This is informative rather than a failure: at full inhibition p53 quickly becomes high enough to saturate both Hill transcription functions, leaving little cell-to-cell variation in the shared input. A real regulatory fork can therefore have almost no snapshot correlation when its driver is clamped or saturated. Conditioning cannot reduce a signal that is already absent. On the real data, conditioning MDM2-CDKN1A on an mRNA p53-activity proxy lowers their correlation only slightly (0.43 to 0.38), so it does not establish a direct edge; the proxy does not measure protein-level p53 activity without error (Kuroki &amp; Pearl, 2014). Partial correlation is nonetheless a standard tool for pruning indirect edges in gene networks (de la Fuente et al., 2004; Sch&auml;fer &amp; Strimmer, 2005).",
    "A fundamental limitation applies: these directional estimators require per-cell time series, whereas scRNA-seq destroys each cell at measurement and yields only population snapshots. Direction on real data must therefore come from population moment dynamics (MBI, Section 4.8) or from RNA velocity, rather than from lagged single-cell statistics."]),
]


def section_results():
    html = ['<h2 id="results">4. Results</h2>',
            '<p class="lead">Each subsection presents one figure with a formal caption and an '
            'interpretation. Sample sizes are stated in every step.</p>']
    def render_block(b):
        if b.get("sub"):
            html.append(f'<p><b>{b["sub"]}</b></p>')
        if b.get("samples"):
            html.append(f'<div class="samples">{b["samples"]}</div>')
        html.append('<figure>' + img(b["fig"]) +
                    f'<figcaption><span class="lab">Figure {b["fign"]}.</span> {b["caption"]}</figcaption></figure>'
                    if b.get("caption") else '<figure>' + img(b["fig"]) + '</figure>')
        for para in b["interp"]:
            html.append(f'<p>{para}</p>')

    for s in RESULTS:
        if s["title"]:
            html.append(f'<h3 id="{s["sid"]}">{s["num"]}&nbsp; {s["title"]}</h3>')
        if "blocks" in s:
            for b in s["blocks"]:
                render_block(b)
        else:
            if s.get("samples"):
                html.append(f'<div class="samples">{s["samples"]}</div>')
            html.append('<figure>' + img(s["fig"]) +
                        f'<figcaption><span class="lab">Figure {s["fign"]}.</span> {s["caption"]}</figcaption></figure>'
                        if s.get("caption") else '<figure>' + img(s["fig"]) + '</figure>')
            for para in s["interp"]:
                html.append(f'<p>{para}</p>')
    return "\n".join(html)


HTML = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>p53-MDM2 moment-based single-cell analysis</title>
<style>{CSS}</style></head>
<body><div class="wrap">
<header>
  <h1>The p53-MDM2 feedback loop and single-cell moment analysis: a simulation and real-data study</h1>
  <p>Built on the higher-order moment framework of Raharinirina et al. (2021), <i>Patterns</i> 2, 100332.</p>
</header>
<main>

<div class="toc">
<b>Contents</b><br>
<a href="#overview">1. Overview</a> &middot;
<a href="#biology">2. The p53-MDM2 loop and Nutlin-3</a> &middot;
<a href="#methods">3. Methods</a> &middot;
<a href="#results">4. Results (4.1-4.10)</a> &middot;
<a href="#takeaways">5. Discussion and conclusions</a> &middot;
<a href="#refs">References</a>
</div>

<h2 id="overview">1. Overview</h2>
<p class="lead">We build a stochastic model of the p53-MDM2 gene circuit, generate synthetic single-cell
RNA-seq from it, and ask a focused question: <b>can the p53&#8594;MDM2 regulatory interaction be
detected&mdash;and correctly directed&mdash;from mRNA measurements?</b> We answer it first on
simulations and then on experimental data, comparing classical dependence measures with the
moment-based inference method of Raharinirina et al. (2021).</p>
<div class="key"><b>Summary of findings.</b> The regulatory signal is real and detectable, but it
resides in higher-order moments and in dynamics rather than in simple correlation; realistic
technical noise (library size) and cell-cycle covariation confound naive correlation; and only
moment-based inference or time-resolved/conditioned statistics recover the directed edge
p53&#8594;MDM2.</div>

<h2 id="biology">2. The p53-MDM2 feedback loop and Nutlin-3</h2>
<p>p53 is a tumour-suppressor transcription factor whose abundance is governed by a two-arm loop
with the oncoprotein MDM2:</p>
<ul>
<li><b>Arm 1 (transcriptional, positive):</b> p53 activates transcription of the <i>MDM2</i> gene,
raising MDM2 mRNA and then MDM2 protein.</li>
<li><b>Arm 2 (post-translational, negative):</b> MDM2, an E3 ubiquitin ligase, targets p53 for
proteasomal degradation, lowering p53.</li>
</ul>
<div class="svgbox">{SCHEMATIC}
<div style="color:#6b7280;font-size:14px;margin-top:6px">p53 turns MDM2 on; MDM2 turns p53 off. Nutlin-3 blocks arm 2.</div></div>
<p>The two arms constitute a negative feedback loop. Because arm 1 passes through an mRNA
intermediate, the loop carries an intrinsic delay, and a delayed negative feedback loop tends to
oscillate; deterministically the oscillation is damped, but molecular noise renders it
self-sustaining&mdash;the p53 pulses observed in single cells (Lahav et al., 2004; Geva-Zatorsky
et al., 2006; Proctor &amp; Gray, 2008).</p>
<p><b>Nutlin-3</b> and its clinical analogue <b>idasanutlin</b> occupy the p53-binding pocket of
MDM2 (Vassilev et al., 2004), abolishing arm 2. p53 accumulates; because arm 1 remains intact,
elevated p53 keeps driving MDM2 transcription, yet the additional MDM2 protein cannot restrain
p53&mdash;the loop is opened. Since the drug acts through p53, it activates p53 targets only in
TP53 wild-type cells, which provides a biological negative control on real data (Section 4.7).</p>

<h2 id="methods">3. Methods</h2>

<h3>3.1 Stochastic model and the Gillespie algorithm</h3>
<p>Molecular species are integer counts evolving as a Markov-jump process; each reaction <i>j</i>
has a propensity a<sub>j</sub>(x). The Gillespie stochastic simulation algorithm draws the waiting
time to the next event and the identity of that event (Gillespie, 1976):</p>
<div class="formula">tau = (1/a<sub>0</sub>) ln(1/r1),   a<sub>0</sub> = sum_j a<sub>j</sub>(x)
choose j*: smallest j with sum_{{k&le;j}} a<sub>k</sub> &gt; r2 a<sub>0</sub>;   x &#8592; x + s<sub>j*</sub>   (r1,r2 ~ U(0,1))</div>
<p>We reuse the exact engine of Raharinirina et al. (2021) (<code>src/ssa.py</code>) and a
numba-accelerated, verified-equivalent implementation (<code>src/fast_ssa.py</code>).</p>

<h3>3.2 The p53-MDM2 model</h3>
<p>Species: p53 mRNA, p53 protein, MDM2 mRNA, MDM2 protein. p53 transcription is constant; MDM2
transcription is activated by p53 (Hill function); p53 degradation by MDM2 is saturable
(Michaelis-Menten). Nutlin efficacy scales the degradation rate by (1 - nutlin).</p>
<div class="formula">MDM2 transcription (arm 1):  k<sub>txn</sub> P^n / (K<sub>d</sub>^n + P^n)        (P = p53)
p53 degradation (arm 2):     (1 - nutlin) k<sub>deg</sub> D P / (K<sub>m</sub> + P)   (D = MDM2)</div>

<h3>3.3 Moments</h3>
<div class="formula">raw moment:  E[A^l1 B^l2](t) = (1/N) sum_n a_n^l1 b_n^l2      (A = p53 mRNA, B = MDM2 mRNA)
covariance:  Cov(A,B) = E[AB] - E[A]E[B]        skewness: E[(A-mu)^3] / sigma^3</div>

<h3>3.4 Dependence measures</h3>
<table>
<tr><th>Measure</th><th>Captures</th><th>Symmetric?</th></tr>
<tr><td>Pearson / Spearman</td><td>linear / monotone correlation</td><td>yes</td></tr>
<tr><td>Mutual information</td><td>any dependence (kNN estimator)</td><td>yes</td></tr>
<tr><td>HSIC, distance correlation</td><td>any (nonlinear) dependence; zero iff independent</td><td>yes</td></tr>
</table>
<p>These answer "are A and B related?" but not "which regulates which?" nor "directly or through a
third factor?".</p>

<h3>3.5 Technical-noise model (Splatter)</h3>
<div class="formula">library size s ~ LogNormal;   counts ~ NB(mean = true * s, dispersion = BCV)
dropout probability:  pi(mu) = 1 / (1 + exp(-k (ln mu - x0)))</div>
<div class="note">The shared per-cell factor s multiplies every gene and therefore induces
<b>spurious</b> gene-gene correlation (Zappia et al., 2017); see Section 4.6.</div>

<h3>3.6 Direction and conditioning</h3>
<div class="formula">Granger (1-step): does A(t) improve prediction of B(t+1) beyond B(t)?
Transfer entropy: TE(A&#8594;B) = I(B<sub>t+1</sub> ; A<sub>t</sub> | B<sub>t</sub>)   (directional MI)
Partial correlation of B and C given A (dependence NOT explained by A):
  residual form:  rB = B - fit(B~A),  rC = C - fit(C~A);   pcorr(B,C|A) = corr(rB, rC)
  closed form:    pcorr(B,C|A) = ( rBC - rBA rCA ) / sqrt( (1 - rBA^2)(1 - rCA^2) )
     where rBC, rBA, rCA are ordinary Pearson correlations</div>
<p>Intuition: regress out the common driver A from both variables and correlate what remains; if B
and C are associated only because both depend on A, the residuals are uncorrelated (pcorr &asymp; 0).
Conditional MI/HSIC are the nonlinear analogues (condition on A instead of removing a linear fit).
Granger causality is due to Granger (1969); transfer entropy to Schreiber (2000); distance
correlation to Sz&eacute;kely et al. (2007) and HSIC to Gretton et al. (2005). The use of conditioning
(partial correlation / conditional independence) to distinguish direct regulation from common-cause
correlation follows the causal-inference framework of Pearl (2009), applied to gene networks by
de la Fuente et al. (2004) and Sch&auml;fer &amp; Strimmer (2005).</p>

<h3>3.7 Moment-based inference (Raharinirina et al., 2021)</h3>
<p>A library of effective mRNA-level reactions is fitted so that the modelled moment dynamics match
the data, closing higher moments with cubic splines and minimizing a non-linear least-squares
residual; a 2x2 directed network is then read off:</p>
<div class="formula">dm/dt = A(theta) m + B(theta) u(t)      (u = higher moments, from data splines)
theta_hat = argmin_{{theta&ge;0}} sum_d || m(t_d|theta) - m_data(t_d) ||^2</div>

<h3>3.8 Permutation z-score</h3>
<p>Because Pearson, MI, HSIC and dCor lie on different scales, their raw values are not comparable.
To test "is the dependence real?" on a common scale we build a null by shuffling one variable
(destroying any relationship), recompute the statistic B times, and standardize:</p>
<div class="formula">z = ( s_obs - mean(s_null) ) / std(s_null)
p = fraction of permutations with |s_perm - mean| &ge; |s_obs - mean|</div>
<p><b>Interpretation.</b> z is the number of standard deviations by which the measured dependence
exceeds chance. z &approx; 0 means nothing beyond chance; z &gt; 2 corresponds roughly to p &lt; 0.05
(a detection); larger z means stronger evidence. Being unitless, z places all measures on the same
footing (Sections 4.6-4.7).</p>

{section_results()}

<h3 id="r-robust">4.10&nbsp; Robustness and reproducibility</h3>
<p>Because the simulator is stochastic, we separate two notions. <b>Reproducibility</b>: the
numba-parallel Gillespie engine is seeded <i>per cell</i>, so a given seed reproduces the ensemble
bit-for-bit (verified: same seed &rarr; identical arrays; a <code>RECOMPUTE</code> flag forces a fresh
run). <b>Robustness</b>: whether the <i>conclusions</i>&mdash;not the exact numbers&mdash;survive changing
the seed. Averages over thousands of cells (Sections 4.1-4.4) are inherently stable; the single-realization
estimates (the permutation z-scores of Section 4.6 and the MBI network of Section 4.8) are the ones that
need checking. For the permutation test we already average over five noise draws; for MBI we re-ran the
inference over several seeds:</p>
<table>
<tr><th>Condition</th><th>A&rarr;B (p53&rarr;MDM2)</th><th>B&rarr;A (MDM2&rarr;p53)</th></tr>
<tr><td>CLOSED loop</td><td>+1.36 &plusmn; 0.60&nbsp; (0.8, 1.0, 2.2)&nbsp; <b>always positive</b></td><td>0.00 &plusmn; 0.00&nbsp; <b>always zero</b></td></tr>
<tr><td>NUTLIN-3</td><td>+1.03 &plusmn; 0.40&nbsp; (1.4, 0.5, 1.2)&nbsp; <b>always positive</b></td><td>&minus;0.26 &plusmn; 0.18&nbsp; (&minus;0.4, 0, &minus;0.3)&nbsp; unstable</td></tr>
</table>
<p>(3 seeds, N = 2500 cells; the main Section 4.8 figure uses N = 8000.) The <b>direction</b> p53&rarr;MDM2
is robust: A&rarr;B is positive in every seed and in both conditions (including Nutlin, where the
transcriptional arm remains intact), while a stable reverse edge never appears (B&rarr;A is exactly zero in
the closed loop and only a small, sign-unstable value under Nutlin&mdash;confirming that the reverse edge
seen in a single Nutlin run is an identifiability artefact, not a real edge). The <i>magnitude</i> of the
edge, by contrast, varies from seed to seed (A&rarr;B ranges 0.5-2.2). The practical lesson: report the
inferred <i>direction</i>, which is trustworthy, rather than the precise edge weight, which is not.</p>

<h2 id="takeaways">5. Discussion and conclusions</h2>
<div class="key">
<ol>
<li>Intrinsic noise converts the deterministically-damped p53-MDM2 loop into sustained p53 pulses;
Nutlin-3 opens the loop, so p53 accumulates and oscillations cease (Sections 4.1-4.2).</li>
<li>The regulatory signal resides in higher-order moments (covariance, skewness), not in the mean of
the p53 transcript (Sections 4.3-4.4), as argued by Raharinirina et al. (2021).</li>
<li>Realistic technical noise (library size) and cell-cycle covariation create spurious
correlations that powerful dependence measures (Pearson, Spearman, HSIC, distance correlation)
all flag as false positives; only the low-power kNN mutual-information estimator stays near zero.
The proper remedy is library-size normalization, not the choice of statistic (Sections 4.6-4.7).</li>
<li>On experimental MIX-seq data the p53 program (MDM2-CDKN1A co-regulation) is detected only in
TP53 wild-type cells under idasanutlin, with a clean isogenic mutant control (Section 4.7).</li>
<li>Symmetric measures report association only; the moment-based inference of Raharinirina et al.
(2021), and time-resolved/conditioned statistics, recover the directed edge p53&#8594;MDM2 and
distinguish regulation from mere correlation (Sections 4.8-4.9).</li>
<li>Obtaining direction from single cells ultimately requires dynamics (moments) or RNA velocity,
because snapshot scRNA-seq has no per-cell time axis.</li>
</ol>
</div>

<h2 id="refs">References</h2>
<div class="refs"><ol>
<li>Raharinirina N.A., Peppert F., von Kleist M., Sch&uuml;tte C., Sunkara V. (2021). Inferring gene
regulatory networks from single-cell RNA-seq temporal snapshot data requires higher-order moments.
<i>Patterns</i> 2, 100332. Code: github.com/vikramsunkara/ScRNAseqMoments.</li>
<li>Proctor C.J., Gray D.A. (2008). Explaining oscillations and variability in the p53-Mdm2 system.
<i>BMC Systems Biology</i> 2:75.</li>
<li>McFarland J.M., Paolella B.R., Warren A., et al. (2020). Multiplexed single-cell transcriptional
response profiling to define cancer vulnerabilities and therapeutic mechanism of action (MIX-seq).
<i>Nature Communications</i> 11, 4296.</li>
<li>Zappia L., Phipson B., Oshlack A. (2017). Splatter: simulation of single-cell RNA sequencing data.
<i>Genome Biology</i> 18:174.</li>
<li>Vassilev L.T., Vu B.T., Graves B., et al. (2004). In vivo activation of the p53 pathway by
small-molecule antagonists of MDM2. <i>Science</i> 303, 844-848.</li>
<li>Lahav G., Rosenfeld N., Sigal A., et al. (2004). Dynamics of the p53-Mdm2 feedback loop in
individual cells. <i>Nature Genetics</i> 36, 147-150.</li>
<li>Geva-Zatorsky N., Rosenfeld N., Itzkovitz S., et al. (2006). Oscillations and variability in the
p53 system. <i>Molecular Systems Biology</i> 2, 2006.0033.</li>
<li>Gillespie D.T. (1976). A general method for numerically simulating the stochastic time evolution
of coupled chemical reactions. <i>Journal of Computational Physics</i> 22, 403-434.</li>
<li>Granger C.W.J. (1969). Investigating causal relations by econometric models and cross-spectral
methods. <i>Econometrica</i> 37, 424-438.</li>
<li>Schreiber T. (2000). Measuring information transfer. <i>Physical Review Letters</i> 85, 461-464.
(transfer entropy)</li>
<li>Sz&eacute;kely G.J., Rizzo M.L., Bakirov N.K. (2007). Measuring and testing dependence by
correlation of distances. <i>Annals of Statistics</i> 35, 2769-2794. (distance correlation)</li>
<li>Gretton A., Bousquet O., Smola A., Sch&ouml;lkopf B. (2005). Measuring statistical dependence with
Hilbert-Schmidt norms. <i>Algorithmic Learning Theory (ALT 2005)</i>, 63-77. (HSIC)</li>
<li>Pearl J. (2009). <i>Causality: Models, Reasoning, and Inference</i>, 2nd ed. Cambridge University
Press. (confounding, common causes, conditioning)</li>
<li>Kuroki M., Pearl J. (2014). Measurement bias and effect restoration in causal inference.
<i>Biometrika</i> 101, 423-437. (proxies of latent confounders; residual confounding)</li>
<li>de la Fuente A., Bing N., Hoeschele I., Mendes P. (2004). Discovery of meaningful associations in
genomic data using partial correlation coefficients. <i>Bioinformatics</i> 20, 3565-3574.</li>
<li>Sch&auml;fer J., Strimmer K. (2005). A shrinkage approach to large-scale covariance matrix estimation
and implications for functional genomics. <i>Statistical Applications in Genetics and Molecular Biology</i>
4, 32. (partial-correlation gene networks)</li>
</ol></div>
</main>
<footer>Generated by <code>report/build_report.py</code>; figures embedded from <code>figures/</code>.
The <code>src/mbi/</code> code is vendored from the repository of Raharinirina et al. (2021); see its ATTRIBUTION.</footer>
</div></body></html>
'''


def main():
    OUT.write_text(HTML, encoding="utf-8")
    print("Wrote", OUT, f"({OUT.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
