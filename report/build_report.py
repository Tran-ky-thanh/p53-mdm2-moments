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
h4{margin-top:25px;font-size:17px;color:#334155;line-height:1.4}
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
           "the experimental data (Section 4.8) used a single near-saturating dose of 2.5 uM (no dose series was available).",
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
   samples="<b>4000 cells per condition</b>, with all four dynamic species initialized at zero; 301 raw snapshots every 4 min from 0 to 1200 min. Panels A-B show the t = 300 min population. The raw moment series are stored in <code>data/cache_step4_zero_initial.npz</code>.",
   caption="Single-cell joint distribution of p53 and MDM2 mRNA with smoothed density contours and the population mean trajectory (orange), for the closed loop (A) and Nutlin-3 (B); mean counts (C), covariance (D) and skewness (E, F) over time. As in Figure 3 of Raharinirina et al. (2021), the orange trajectory begins at zero expression. Curves are smoothed only for display (Savitzky-Golay; PCHIP for the orange path); all reported statistics and cached values remain the unsmoothed 4000-cell estimates.",
   interp=[
    "In the closed loop the joint cloud is tilted along a positive diagonal (covariance +10.8, correlation +0.43): cells rich in p53 mRNA also tend to be rich in MDM2 mRNA, the fingerprint of an intact transcriptional coupling. Under Nutlin the cloud becomes round and shifts upward (covariance -0.4, correlation -0.01, both effectively zero), with MDM2 mRNA elevated&mdash;dependence is lost even though the marginal level of MDM2 mRNA is higher.",
    "Why did the previous version not begin at zero? It inherited the model's basal starting state (15 TP53 mRNAs, 30 p53 proteins, 20 MDM2 mRNAs and 50 MDM2 proteins). Moreover, 15 is exactly the mean steady-state TP53-mRNA count implied by constitutive transcription and decay: dE[m<sub>P</sub>]/dt = 3.0 - 0.20E[m<sub>P</sub>], which equals zero at E[m<sub>P</sub>] = 15. The old green marker therefore correctly appeared at (15,20), and TP53 mRNA had no initial mean change, but that did not reproduce the visual experiment in Patterns2021. The revised Figure 4 instead starts (TP53 mRNA, p53 protein, MDM2 mRNA, MDM2 protein) = (0,0,0,0). Its expected TP53-mRNA start-up is E[m<sub>P</sub>](t) = 15(1 - exp[-0.20t]), so it rises immediately from zero toward 15 counts. This zero start is a deliberate start-up experiment for Figure 4, not a claim that a treated biological cell contains no p53 or MDM2 at baseline.",
    "Panel C shows that after the start-up transient the two solid p53-mRNA curves coincide, underscoring that the discriminating information lies in the second- and third-order statistics (D-F), not in the first-order mean of p53 mRNA."]),
 dict(sid="r6b", num="4.5", fig="fig_step6b_normalization.png", fign=5,
   title="Realistic technical noise, the library-size confound, and whether preprocessing removes it",
   samples="<b>3000 cells per condition</b> plus <b>300 background genes</b> (a stand-in for the rest of the transcriptome, independent of p53/MDM2 but sharing each cell's library-size factor); snapshot at t = 600 min; Splatter noise with dropout midpoint +1 (about 14% zeros in the two genes of interest). Panels A-D average <b>5</b> independent noise draws (mean &plusmn; SD); panels E-F average <b>3</b> noise draws of a 150-permutation test on a 1200-cell subsample. Cached in <code>data/cache_step6b.npz</code>; script <code>analysis/step6b_normalization.py</code>.",
   caption="Five preprocessing pipelines evaluated against the <b>clean-count ground truth</b> (dotted horizontal lines &mdash; defined in the first paragraph below). (A-D) Pearson correlation, kNN mutual information, normalized HSIC and distance correlation for the closed loop (blue) and Nutlin-3 (red). (E) Permutation z-score in the closed loop, where the dependence is real and a <i>high</i> value is correct. (F) Permutation z-score under Nutlin-3, where the two transcripts are genuinely independent and any value above the z = 2 line is a <i>false positive</i>. Pipelines: raw counts (no preprocessing), log1p without normalization, counts-per-10,000 followed by log1p (the Seurat/Scanpy standard), the same preceded by a quality-control filter that removes the lowest 5% of cells by sequencing depth and by genes detected, and finally CP10k + log1p with log sequencing depth regressed out of each gene.",
   interp=[
    "This section replaces the idealized dropout of the earlier steps with the published <b>Splatter</b> technical-noise model (Zappia et al., 2017; Method 3.5): a per-cell log-normal library-size factor that multiplies every gene, negative-binomial overdispersion, and expression-dependent logistic dropout. The shared library-size factor is the key troublemaker&mdash;because it scales all genes in a cell together, it induces a <i>spurious</i> positive correlation between transcripts that are in fact independent. This is the well-known library-size confound and a concrete example of the pre-processing pitfalls emphasized by Raharinirina et al. (2021). The leftmost <b>raw-counts</b> column of every panel shows it directly: under Nutlin, where p53 and MDM2 are genuinely independent, the correlation-type measures still report non-zero values (dCor 0.098, panel D) and, worse, the permutation test calls them significant (panels E-F). The rest of this section asks whether standard preprocessing removes the confound&mdash;and at what cost to the real signal.",
    "The <b>raw-counts</b> pipeline feeds noisy counts straight into the four statistics, which no real analysis pipeline would do. Two design choices make the comparison a fair test, and each answers a natural question about the figure. <b>What are the dotted lines?</b> They are the <i>clean-count ground truth</i>: each statistic computed on the simulated SSA counts <i>before</i> any technical noise is added&mdash;the true dependence between p53 and MDM2 mRNA that a perfect pipeline would recover. A bar sitting on its dotted line means that pipeline recovered the truth; a bar far from it&mdash;especially a tall bar under Nutlin, where the truth is near zero&mdash;is reporting an artefact. <b>Why add 300 genes?</b> Library-size normalization divides each cell's counts by that cell's <i>total</i> transcript count, so it is meaningless with only two genes: dividing two numbers by their own sum is circular. Our simulation produces only p53 and MDM2 mRNA, so we add 300 &lsquo;spectator&rsquo; genes that stand in for the rest of the transcriptome (a real cell expresses ~20,000). They are statistically independent of p53 and MDM2, but Splatter multiplies them by the <i>same</i> per-cell library-size factor as the two genes of interest&mdash;so they reproduce the confounder exactly as it arises in real data and give normalization a realistic total to divide by. We then apply the same steps used by Seurat (Hao et al., 2021) and Scanpy (Wolf et al., 2018).",
    "<b>How to read the permutation z-score (panels E-F).</b> The <i>value</i> of a measure (say dCor = 0.05) and its <i>permutation z-score</i> answer two different questions, and it is worth separating them before reading the results. The value asks <b>how strong</b> the association is; the z-score asks <b>whether there is any association at all</b>, by comparing the observed statistic against a null distribution built from randomly shuffling the two genes against each other (Method 3.8). The z-score is a detector, not a measure of how trustworthy the number is: z &gt; 2 means &lsquo;a relationship was detected&rsquo;, z &lt; 2 means &lsquo;none detected&rsquo;. This is why, under Nutlin, a <i>small</i> z is the <b>correct</b> result&mdash;the two genes really are independent, so the detector should stay quiet, exactly as a smoke alarm should stay silent in a room with no fire. Running the permutation test on the open loop is precisely the false-alarm check: we <i>want</i> z &lt; 2 there, and obtaining it (after normalization, panel F) is what demonstrates that the pipeline no longer invents associations.",
    "<b>A large z does not mean a strong relationship.</b> Because the z-score grows with sample size even for a tiny correlation, it measures <i>confidence</i>, not <i>strength</i>. The raw-count Nutlin bars make this concrete: dCor is only 0.098 (a small value) yet z = +3.0 (a confident false positive), and with more cells that same r &asymp; 0.09 artefact climbs to z = +9 (Figure 6). A small value therefore does <i>not</i> by itself protect against a false positive&mdash;only the z-score, read against a proper negative control, separates &lsquo;small and real&rsquo; from &lsquo;small and spurious&rsquo;. The two numbers are complementary and the report always shows both: the value (panels A-D) is the effect size, the z-score (panels E-F) says whether to believe it.",
    "<b>Normalization removes the artefact, and it does so completely.</b> Panel F is the decisive one. On raw counts the Nutlin-3 condition&mdash;where the truth is independence&mdash;yields significant false positives for distance correlation (z = +3.0), Pearson (z = +2.0) and HSIC (z = +2.1). After CP10k + log1p every one of them collapses to the null band (dCor z = +0.5, Pearson z = -0.9, HSIC z = -0.3). Panels A and D show the same result on the effect-size scale: raw counts report dCor = 0.098 under Nutlin where the clean-count truth is 0.046, whereas CP10k + log1p returns 0.045, matching the truth to within the noise. This directly confirms that normalization is the remedy for the library-size confound (the raw-counts column) and, importantly, it is the <i>normalization</i> that matters: log1p alone (second pipeline) only halves the artefact, and the additional QC filter and depth regression change little beyond what CP10k already achieved.",
    "<b>But the honest accounting is less comfortable.</b> Panel E shows that normalization also removes most of the apparent detection in the closed loop, where the dependence is real: dCor falls from z = +4.9 to z = +0.1 and Pearson from +3.7 to +0.9. The correct interpretation is not that normalization destroyed a real signal, but that the raw-count z-scores at the left of panels E-F were largely measuring the library-size artefact in the first place. Comparing effect sizes makes this explicit: of the closed-loop raw dCor of 0.121, the Nutlin condition&mdash;which contains artefact only&mdash;already accounts for 0.098. Once that shared component is removed, the genuine residual is small (dCor 0.049 against a clean-count truth of 0.368), because Splatter's negative-binomial overdispersion (BCV = 0.5) and Poisson sampling independently corrupt each gene and attenuate the true correlation. A dropout sweep confirms that this attenuation is driven mainly by the overdispersion rather than by dropout: even at 1% zeros the normalized closed-loop dCor is only 0.083. The practical consequence is that with 3000 cells this residual signal no longer reaches significance.",
    "<b>Can more cells recover it?</b> Figure 6 answers this with a cell-number sweep, and the answer is encouraging for the normalized pipeline and alarming for the unnormalized one. After CP10k + log1p the closed-loop correlation is genuine but small (r between 0.016 and 0.033 across sample sizes), so its significance accumulates slowly, reaching z = +3.6 only at 12,000 cells; the crucial point is that its negative control moves in the opposite direction over the same range (Nutlin z = +0.2 down to -1.0). Normalization therefore does not destroy the real signal&mdash;it makes the signal expensive but honest, and around 10<sup>4</sup> cells it becomes detectable with the correct null. This also settles the interpretation of panel E in Figure 5: the collapse of the z-scores there reflects the removal of the artefact plus the modest sample size of 3000 cells, not the loss of a real effect.",
    "The raw-count curves carry the warning. Under Nutlin, where the correct answer is no association at all, significance climbs from z = +3.9 at 1500 cells to z = +9.5 at 12,000 while the underlying correlation stays essentially constant at r &asymp; 0.09. Because the library-size factor induces a real, reproducible population-level correlation&mdash;a bias rather than sampling noise&mdash;a larger sample estimates it ever more precisely and reports it as ever more significant. Stated plainly: <b>on unnormalized counts a bigger experiment does not protect against this artefact, it entrenches it.</b> Note also that the raw closed-loop and raw Nutlin curves are separated by only a constant offset, which is exactly what one expects when both are dominated by the same shared factor.",
    "<b>An unexpected result for mutual information.</b> In the raw-count analyses above the kNN (Kraskov et al., 2004) mutual-information estimator was the weakest detector, and we attributed this to low power on discrete count data. Panels B and E refine that explanation. After CP10k + log1p with log depth regressed out, the estimator recovers MI = 0.070 nats in the closed loop against a clean-count truth of 0.087, and it becomes the <i>only</i> measure that still exceeds z = 2 there (z = +3.0) while remaining at the null under Nutlin (z = +0.5, panel F). The likely mechanism is that raw UMI counts are heavily tied&mdash;many cells share the identical integer value&mdash;and the Kraskov estimator's nearest-neighbour distances degenerate under ties. Residualization maps the counts to continuous values, breaks the ties, and restores the estimator. We report this as a single-configuration observation rather than a general recommendation: the closed-loop and Nutlin values are separated by more than two standard deviations across the five noise draws (0.070 &plusmn; 0.018 versus 0.022 &plusmn; 0.008), but confirming it would require a systematic benchmark across noise levels and cell numbers.",
    "Note that <i>scaling</i> in the usual sense&mdash;z-scoring each gene&mdash;cannot appear as a separate pipeline here, because it is an affine transformation and Pearson, Spearman and distance correlation are all invariant to it; normalized HSIC with the median-distance bandwidth heuristic is invariant as well. Scaling matters for downstream PCA and clustering, not for these pairwise dependence statistics."]),
 dict(sid="r6c", num=None, fig="fig_step6b_ncells.png", fign=6,
   title=None,
   samples="One <b>12,000-cell</b> ensemble per condition, subsampled without replacement to 1500, 3000, 6000 and 12,000 cells so that only the sample size differs; each point averages <b>3</b> independent Splatter noise draws with the same 300 background genes as Figure 5.",
   caption="Significance of the Pearson correlation between p53 and MDM2 mRNA as a function of the number of cells, on raw counts (solid, filled markers) versus CP10k + log1p (dashed, open markers), for the closed loop (blue) and Nutlin-3 (red). Significance is the Fisher transform z = arctanh(r)&radic;(n-3); the dotted line marks z = 2 (approximately p &lt; 0.05). Under Nutlin-3 the two transcripts are independent by construction, so the red solid curve measures nothing but the library-size artefact.",
   interp=[
    "The red solid curve is the most important one. It tracks a condition in which the correct answer is <i>no association</i>, yet its significance rises steadily with sample size, from z = +3.9 at 1500 cells to z = +9.5 at 12,000, while the correlation it is testing barely moves (r &asymp; 0.10 falling to 0.09). The shared per-cell library-size factor induces a real, reproducible correlation in the observed counts; it is a bias, not sampling noise, so collecting more cells estimates it more precisely instead of averaging it away. Significance answers only &lsquo;is this correlation distinguishable from zero?&rsquo;, and for an artefact of this kind the honest answer is yes, increasingly so. The blue solid curve sits above it but rises in parallel, consistent with both being dominated by the same shared factor.",
    "The dashed curves show what normalization buys. The Nutlin control (red, dashed) stays in the null band throughout and in fact drifts slightly negative (z = +0.2 at 1500 cells to -1.0 at 12,000), confirming that the artefact is gone rather than merely reduced. The closed loop (blue, dashed) starts indistinguishable from it but separates as cells accumulate, crossing the z = 2 line between 6000 and 12,000 cells and reaching z = +3.6 at 12,000. The residual true effect is small, r between 0.016 and 0.033, and because only three noise draws contribute to each point these individual values carry appreciable uncertainty; the reliable reading is the divergence between the two dashed curves, not any single z-value.",
    "Two practical conclusions follow. First, library-size normalization must be applied <i>before</i> any pairwise dependence statistic, and the size of the experiment is no substitute for it&mdash;an unnormalized analysis of 12,000 cells is more confidently wrong than one of 1500. Second, because the surviving effect sizes are small, significance alone is not a sufficient criterion; a negative control condition&mdash;Nutlin-3 here, and the TP53-mutant cell lines in the experimental data of Section 4.6&mdash;is what distinguishes a weak real signal from a residual artefact. The alternative to spending 10<sup>4</sup> cells on a single snapshot correlation is to aggregate information over time instead of over cells, which is the moment-dynamics argument of Raharinirina et al. (2021) taken up in Section 4.7."]),
 dict(sid="r8", num="4.6", title="Validation on experimental scRNA-seq: MIX-seq DMSO versus idasanutlin",
   blocks=[
    dict(sub="Why not the direct (p53, MDM2) pair, and why does the sign look reversed versus the simulation?",
      fig="fig_step7d_direct_pair.png", fign=7,
      samples="WT line LNCaP, 6 h (n=152 treated). dCor of TP53-MDM2 vs MDM2-CDKN1A.",
      caption="For the wild-type line: (A) mean expression of TP53, MDM2 and CDKN1A in DMSO vs idasanutlin; "
              "(B) scatter of the direct loop pair TP53 mRNA vs MDM2 mRNA; (C) dCor of TP53-MDM2 vs the "
              "p53-target pair MDM2-CDKN1A in each condition.",
      interp=[
       "Our ideal observable is the loop pair (p53, MDM2), and the simulation uses exactly that. On real scRNA-seq, however, <i>TP53 mRNA</i> is a poor readout of p53 <i>activity</i>: p53 is controlled post-translationally, so its transcript barely moves and does not track the protein that actually drives MDM2 transcription. Panel (A) confirms this directly&mdash;TP53 mRNA barely changes with the drug (it even dips slightly), while MDM2 and CDKN1A are strongly induced. Panels (B)-(C) show the direct TP53-MDM2 scatter has little structure (dCor 0.07 in DMSO, 0.10 under idasanutlin) and stays weak in both conditions, whereas the p53-<i>target</i> pair MDM2-CDKN1A jumps from 0.12 to 0.38. On mRNA-only data the loop's regulatory coherence is therefore best read from p53-target co-expression, not from the nominal (p53, MDM2) pair&mdash;an intrinsic limitation of measuring a post-translationally regulated hub at the mRNA level. We use MDM2-CDKN1A for the rest of this section.",
       "The <i>direction</i> of the drug effect on that pair also differs from the simulation, for a principled reason: cell-to-cell covariance requires both an active regulatory link <b>and</b> heterogeneity in the upstream driver. In the simulation the baseline (no drug) is an actively cycling loop, so covariance is high at baseline and Nutlin&mdash;by clamping p53&mdash;removes it. In the experiment the baseline (DMSO) is quiescent (low, uniform p53 activity), so target co-expression is low, and it is the drug that switches p53 on heterogeneously, creating the co-expression. Same principle, opposite baseline."]),
    dict(sub=None, fig="fig_step7_realdata.png", fign=8,
      samples="Experimental data from MIX-seq (McFarland et al., 2020): <b>1886 (DMSO) and 2634 (idasanutlin) QC-passing singlets</b> across a pool of 24 cell lines; the TP53 wild-type line LNCaP contributes n &approx; 100 / 152 cells; permutation tests use B = 200-300.",
      caption="Application to real scRNA-seq: DMSO (closed loop) versus idasanutlin, a clinical MDM2 inhibitor of the Nutlin family (open loop), after CP10k + log normalization. (A) induction of p53 target genes in a TP53-wild-type line but not a mutant line; (B, C) MDM2-CDKN1A joint distributions in the wild-type and mutant lines; (D) all five measures on the wild-type MDM2-CDKN1A pair; (E) the 2x2 wild-type/mutant biological control; (F) change in dCor upon treatment across gene pairs.",
      interp=[
       "Idasanutlin sharply induces canonical p53 targets (MDM2, CDKN1A, FDXR, ...) in the TP53-wild-type line LNCaP but not in the TP53-mutant line, as expected for a p53-dependent drug (McFarland et al., 2020). Concomitantly, MDM2 and CDKN1A&mdash;two direct p53 targets&mdash;become co-expressed under treatment: dCor rises from 0.12 (not significant) to 0.38 (permutation z &approx; 9). The 2x2 control (E) isolates the p53-dependent effect: only the wild-type + idasanutlin condition shows strong co-regulation, while the drug-treated mutant line does not.",
       "A caveat is made explicit in panel F: among all gene pairs, only MDM2-CDKN1A gains dependence upon treatment, whereas housekeeping pairs actually lose covariation&mdash;idasanutlin arrests the cell cycle in wild-type cells, removing a global source of gene-gene covariance. Housekeeping pairs are therefore <i>not</i> a perfect negative control in single-cell data; the cleanest control is the isogenic comparison against a non-responding (mutant) line."]),
    dict(sub="Does the effect generalize across cell lines?",
      fig="fig_step7e_across_lines.png", fign=9,
      samples="All <b>22 MIX-seq cell lines with &ge; 40 QC-passing cells in both DMSO and idasanutlin</b> (6 h). "
              "For each line: p53 response = mean induction of MDM2 + CDKN1A (idasanutlin minus DMSO, log-norm); "
              "dCor gain = distance correlation of (MDM2, CDKN1A) in idasanutlin minus in DMSO. Marker area is "
              "linearly proportional to the total number of QC-passing cells for that line across both conditions.",
      caption="One point per cell line: p53 response (x-axis) versus the gain in MDM2-CDKN1A co-expression "
              "under idasanutlin (y-axis), with a linear trend fitted across all 22 lines. Point area represents "
              "the number of cells retained after preprocessing in DMSO plus idasanutlin; the size legend gives "
              "the minimum, median and maximum totals.",
      interp=[
       "There is a positive association (r = +0.44): lines whose p53 targets are induced more strongly by idasanutlin also tend to gain more MDM2-CDKN1A co-expression, consistent with the mechanism above. LNCaP, the wild-type line used throughout this section, sits at the extreme of both axes and is the cleanest example in the panel; CCFSTTG1, NCIH226 and DKMG (moderate p53 response) also show a modest positive gain.",
       "The relationship is noisy, however, and we report it as such rather than overstating it. Two lines with almost no p53 response (RCM1, BT549) show large <i>negative</i> dCor changes, most likely small-sample or line-specific noise rather than a p53-related effect, and they visibly pull down the correlation. With only 22 lines and per-line cell counts of a few hundred, this trend should be read as supportive but not as strong independent confirmation; it motivates, rather than replaces, the single-line analysis with its 2x2 biological control."]),
    dict(sub="Timepoint dependence: 6 h versus 24 h of idasanutlin",
      fig="fig_step7c_timepoints.png", fign=10,
      samples="Idasanutlin was applied at a single saturating concentration (<b>2.5 uM</b>) and assayed at "
              "<b>6 h and 24 h</b> (MIX-seq has no idasanutlin dose series). WT line LNCaP: 152 vs 102 treated cells.",
      caption="For the TP53 wild-type line, (A) induction of p53 target genes (idasanutlin minus DMSO) at "
              "6 h versus 24 h, and (B) MDM2-CDKN1A co-expression (dCor) at the two timepoints.",
      interp=[
       "The two readouts move in <i>opposite</i> directions with time. The <b>mean induction</b> of p53 targets is generally <b>larger at 24 h</b> (summed log-fold induction 6.6 at 6 h vs 9.4 at 24 h), reflecting continued accumulation of target mRNA. In contrast, the <b>cell-to-cell co-expression</b> of MDM2 and CDKN1A is <b>stronger at 6 h</b> (dCor 0.38) than at 24 h (0.26).",
       "This is expected: at 6 h the response is still in its heterogeneous transient, so p53 activity&mdash;and hence its targets&mdash;varies markedly from cell to cell, maximizing the covariation our measures detect; by 24 h the response has largely saturated and homogenized (and fewer cells are recovered, consistent with arrest/apoptosis), so the mean is higher but the covariation is lower. The rest of this section uses the 6 h timepoint, which gives the cleaner co-regulation signal. Note that 'efficacy' in our simulation (Section 4.2) is a mechanistic fraction-of-binding-blocked in [0,1], whereas 2.5 uM is a single, near-saturating <i>concentration</i>; the two are related monotonically but are not the same quantity."]),
   ]),
 dict(sid="r9", num="4.7", fig="fig_step8_mbi.png", fign=11,
   title="Recovering directed regulation from moments: non-linear moment-based inference (Raharinirina et al., 2021)",
   samples="Panels A-B use one <b>8000-cell simulation per condition</b> to show the fitted mean moment dynamics. Panels C-E summarize <b>5 independent seeds &times; 2500 cells per condition</b>; each replicate uses 25 moment snapshots from 20 to 980 min after discarding the first 20 min and subsampling every 40 min. Raw moments up to order 3-4 are used. "
           "<b>Note on noise:</b> the input is <i>clean</i> SSA counts&mdash;they carry the <i>intrinsic</i> (molecular) noise that MBI actually exploits, but <i>no</i> technical scRNA-seq noise (dropout, library size) is applied here. This is a best-case test, as in Raharinirina et al. (2021).",
   caption="Non-linear moment-based inference (MBI) of Raharinirina et al. (2021), applied to simulated moment time-courses of the (p53 mRNA, MDM2 mRNA) pair (original code, vendored in src/mbi). A-B: one 8000-cell example fit to the mean dynamics. C-D: inferred 2x2 regulatory networks, shown as mean &plusmn; SD across five independent 2500-cell simulations. E: the two directed edge strengths, again shown as mean &plusmn; SD across seeds.",
   interp=[
    "Unlike the symmetric dependence measures of the previous steps, MBI returns a <i>directed</i> network. In the closed loop it consistently recovers the transcriptional edge p53 &#8594; MDM2: across five independent simulations A&#8594;B = +1.66 &plusmn; 0.75 and is positive in every seed (0.84, 1.08, 2.76, 1.82, 1.78). The reverse mRNA-level edge is absent: B&#8594;A = +0.01 &plusmn; 0.02. This is biologically expected, because the negative arm operates post-translationally on p53 protein and is invisible as an MDM2-mRNA &#8594; TP53-mRNA interaction.",
    "Under Nutlin-3 MBI still recovers a positive p53 &#8594; MDM2 edge (A&#8594;B = +0.99 &plusmn; 0.60; positive in all five seeds). This is the correct biology: Nutlin blocks the MDM2&ndash;p53 <i>protein</i> interaction, not the transcriptional arm, so p53 continues to drive MDM2 transcription. The reverse estimate is weaker and less stable (B&#8594;A = &minus;0.27 &plusmn; 0.22). With the loop open, p53 is high and near-static and MDM2 mRNA sits on a plateau, so the moment time-courses contain less dynamic information for the fit. MBI is therefore most reliable in the informative, oscillating closed-loop regime; a near-static system gives noisier edge weights.",
    "The practical robustness result is that the <i>direction</i> p53 &#8594; MDM2 survives changing random seeds, while the exact edge magnitude does not. This is why the report interprets MBI qualitatively as directional support rather than as a calibrated biochemical rate. Reproducibility is handled separately: the numba-parallel Gillespie engine is seeded per cell, so the same seed reproduces the ensemble exactly, and <code>RECOMPUTE=1</code> intentionally forces a fresh simulation.",
    "Finally, a scope note: MBI works on moment time-courses and implicitly assumes the <i>observed</i> moments equal the <i>true</i> moments. It de-noises only by averaging over many cells and by spline-smoothing higher moments, which suppresses sampling noise&mdash;but it does not model technical scRNA-seq artefacts. Because dropout and library-size shifts distort exactly the higher-order moments MBI relies on (Section 4.5), such technical noise would bias the inference; the clean, intrinsic-noise-only setting here is deliberately a best case."]),
 dict(sid="r10", num="4.8", fig="fig_step9_directional.png", fign=12,
   title="Endowing dependence measures with direction and conditional specificity",
   samples="<b>4000 cells per condition</b> for the directional analysis (fluctuation window t &ge; 200 min). Panels D-E use a separate <b>4000-cell, six-species p53-MDM2-CDKN1A Gillespie simulation</b>. Each independent cell starts from (TP53 mRNA, p53 protein, MDM2 mRNA, MDM2 protein, CDKN1A mRNA, CDKN1A protein) = (15, 30, 20, 50, 15, 40) molecules and is recorded at <b>21 snapshots</b>, every 4 min from 0 to 80 min, with Nutlin efficacy = 1. The clean simulation contains intrinsic reaction noise only. Panel E additionally applies Splatter-style technical noise to the observed MDM2 and CDKN1A transcript counts: log-normal library-size variation, negative-binomial overdispersion and expression-dependent dropout. The realized zero fractions after this noise are 5.7% for MDM2 and 5.7% for CDKN1A. The random seed is 93; the asynchronous-snapshot seed is 101; the Splatter seed is 202. The complete output is saved in <code>data/cache_step9_three_gene_nutlin.npz</code>; all 20 kinetic parameters are tabulated in Method 3.2. Real-data partial correlation uses LNCaP (n &approx; 152).",
   caption="(A) lagged cross-correlation and (B) lagged distance correlation of the detrended p53/MDM2 mRNA "
           "fluctuations, shown for tau &ge; 0 as two directional curves; at positive lag the p53&rarr;MDM2 curve "
           "dominates and peaks near tau &asymp; +4 min while MDM2&rarr;p53 decays, i.e. p53 fluctuations precede "
           "MDM2. (C) one-step Granger causality in both directions; (D) mean TP53-mRNA, p53-protein, MDM2-mRNA and CDKN1A-mRNA "
           "response trajectories from an explicit six-species simulation under full Nutlin; (E) marginal-versus-partial "
           "correlations from the clean synthetic asynchronous Nutlin snapshot, the Splatter-noised version of the same "
           "snapshot, and real idasanutlin data.",
   interp=[
    "Two ingredients convert the symmetric measures of Sections 4.5-4.6 into causal ones. <b>Direction requires time</b>: using per-cell trajectories, the lagged cross-correlation and lagged distance correlation peak at positive lag for p53&#8594;MDM2 (A leads B). In panel C, one-step Granger causality is also strongly asymmetric in the dynamic closed loop: p53&#8594;MDM2 = 0.206699, whereas MDM2&#8594;p53 = 2.23 &times; 10<sup>&minus;6</sup>. Under Nutlin the corresponding values are 9.75 &times; 10<sup>&minus;8</sup> and 1.16 &times; 10<sup>&minus;6</sup>. Thus the purple bars are not mathematically zero; the earlier labels displayed 0.000 only because they were rounded to three decimals. They are nevertheless negligible on the scale of the recovered forward closed-loop signal. The open-loop mRNA trajectories rapidly settle onto a nearly static plateau, leaving little predictive fluctuation. (Because the lagged measures are symmetric, C<sub>B&#8594;A</sub>(&tau;)=C<sub>A&#8594;B</sub>(&minus;&tau;), it is enough to show <i>positive</i> lags only: over &tau; &ge; 0 the two directional curves are genuinely distinct, and the p53&#8594;MDM2 curve exceeding MDM2&#8594;p53 there means p53 fluctuations <i>precede</i> MDM2&mdash;the temporal-precedence signature of p53&#8594;MDM2. The ~4-min peak reflects the transcription/translation delay.) A further advantage of the <i>signed</i> cross-correlation (panel A) is the <b>negative lobe near &tau; &asymp; 16 min</b>: after p53 drives MDM2 up, the negative-feedback arm pulls the system back, giving an anti-correlation (overshoot) at roughly half the loop's oscillation timescale. This lobe is <i>consistent with</i> the closed negative-feedback loop&mdash;it is absent under Nutlin (open loop; grey dashed curve) and, being a change of <i>sign</i>, cannot appear in the non-negative distance correlation of panel B.",
    "<b>Panels D and E use an explicit p53-MDM2-CDKN1A simulation.</b> It contains mRNA and protein for all three genes and 13 reaction channels, simulated exactly with the Gillespie direct method. All 4000 cells use the same initial molecular counts and kinetic constants but have independent reaction-event histories, so their variation is intrinsic stochastic variation. Here <i>open loop</i> means that Nutlin removes only the reverse inhibitory arm, MDM2 protein &#8866; p53; it does not inhibit p53 as a transcription factor. At efficacy 1, the effective MDM2-dependent p53 degradation constant is exactly zero, while the forward branches p53 &#8594; MDM2 mRNA and p53 &#8594; CDKN1A mRNA remain. This gives the common-driver fork MDM2 mRNA &#8592; p53 protein &#8594; CDKN1A mRNA. Panel D plots log<sub>2</sub>[(mean count at t + 1)/(mean count at t = 0 + 1)]. TP53 mRNA remains close to baseline because its transcription is constitutive in the model, whereas p53 protein rises continuously to more than 6 log2-fold above baseline after MDM2-mediated degradation is blocked. MDM2 and CDKN1A mRNAs also increase and then plateau because both forward transcriptional branches remain active. This separation between a nearly unchanged TP53 transcript and strongly accumulating p53 protein illustrates why TP53 mRNA is a poor proxy for p53 activity in scRNA-seq. The rising MDM2 transcript cannot restore feedback because Nutlin prevents its protein product from degrading p53. CDKN1A is a canonical p53 target associated with cell-cycle arrest; this minimal model simulates its expression but does not simulate the cell cycle itself.",
    "Panel E converts the trajectories into a synthetic asynchronous scRNA-seq-like snapshot: each of the 4000 cells contributes exactly one observation, sampled uniformly from the 20 post-treatment times (4-80 min). In the clean snapshot, response-phase heterogeneity makes MDM2 and CDKN1A strongly correlated (r = 0.602). The conditioning covariates are two target-specific, exponentially weighted histories of p53-protein Hill activity, using each target's Hill constant, Hill coefficient and mRNA-decay rate. They represent the shared p53 input accumulated before the sampled time. Regressing both target mRNAs on these covariates and correlating the residuals reduces the partial correlation to -0.001. Thus the clean tall marginal bar reflects the common p53 fork, rather than a direct MDM2-CDKN1A edge.",
    "After adding Splatter-style technical noise to the two observed target transcripts, the marginal correlation falls to r = 0.379 but the partial correlation remains high (r = 0.338). This happens because the added library-size factor and dropout affect the observed MDM2 and CDKN1A counts directly, whereas the conditioning variable is the biological p53-protein exposure. Conditioning removes the biological common driver but cannot remove a technical common driver that is not included in the model. The real idasanutlin data now sit beside this noisy synthetic control in the same panel: conditioning MDM2-CDKN1A on an mRNA p53-activity proxy lowers their correlation only slightly (0.426 to 0.376). The similarity between the noisy synthetic bars and the real-data bars supports the interpretation that real partial correlations can remain high because of imperfect p53-activity measurement plus residual technical and biological confounding, rather than because MDM2 directly regulates CDKN1A. Partial correlation is nonetheless a standard tool for pruning indirect edges in gene networks (de la Fuente et al., 2004; Sch&auml;fer &amp; Strimmer, 2005). Figure 13 tests this interpretation directly by asking whether library-size <i>normalization</i>&mdash;which removes the technical common driver&mdash;makes conditioning work again.",
    "A fundamental limitation applies: these directional estimators require per-cell time series, whereas scRNA-seq destroys each cell at measurement and yields only population snapshots. Direction on real data must therefore come from population moment dynamics (MBI, Section 4.7) or from RNA velocity, rather than from lagged single-cell statistics."]),
 dict(sid="r10b", num=None, fig="fig_step9_normalization.png", fign=13,
   title=None,
   samples="Follow-up on panel E of Figure 12. Synthetic: the same 4000-cell asynchronous Nutlin snapshot, re-noised with Splatter plus <b>300 background genes</b> so a library size exists (seed 202). Real: the LNCaP idasanutlin cells (n = 152). Each is evaluated at three preprocessing levels; the synthetic conditions on the biological p53-protein exposure, the real data on an mRNA p53-target proxy normalized the same way as the targets. Cached in <code>data/cache_step9_norm.npz</code>.",
   caption="Does library-size normalization remove the residual partial correlation of Figure 12E? "
           "(A) Splatter-noised synthetic snapshot and (B) real LNCaP idasanutlin cells, each on raw counts, "
           "log1p only, and CP10k + log1p. Blue = marginal correlation of MDM2 and CDKN1A; red = partial "
           "correlation after conditioning on the p53 driver. Dotted lines in (A) mark the clean-count ideal "
           "(marginal &asymp; 0.60, partial &asymp; 0).",
   interp=[
    "In the synthetic system the truth is known, so it calibrates the test. Conditioning on the real p53 driver removes the MDM2-CDKN1A correlation entirely on clean counts (partial &minus;0.001). On raw Splatter counts conditioning <i>fails</i> (partial +0.354), because the shared library-size factor is a technical common driver absent from the conditioning set. Normalization fixes exactly this: <b>log1p alone</b> only roughly halves the residual (partial +0.180), whereas <b>CP10k + log1p</b> brings it back to +0.031&mdash;essentially the clean-count ideal. This confirms the panel-E interpretation directly: the residual there was the technical library-size confound, and removing it lets conditioning recover the truth. (The raw value here, 0.354, differs trivially from the 0.338 quoted for Figure 12E only because this analysis adds 300 background genes so that a library size can be defined.)",
    "The real data behave differently, and that is the key result. <b>log1p alone changes nothing</b> (partial 0.439 &rarr; 0.440), and even full <b>CP10k + log1p</b> lowers the partial correlation only modestly, to 0.376&mdash;nowhere near the synthetic's 0.031. So on real idasanutlin data the residual partial correlation is <i>not</i> a library-size artefact: the very normalization that erased it in the synthetic control (A) leaves it almost untouched here (B). What remains must therefore be biological or measurement-driven&mdash;an imperfect mRNA proxy for p53 activity, plus genuine co-regulation or shared confounders beyond the single p53 input we condition on. The comparison sharpens the earlier conclusion: normalization is necessary and sufficient to remove the <i>technical</i> confound, but a stubbornly high real-data partial correlation reflects the limits of inferring regulation from mRNA snapshots, not a preprocessing failure. It also answers the practical question &lsquo;does log1p help?&rsquo; cleanly: log1p by itself is not enough (it does not remove the shared library-size factor); the counts-per-10,000 step is what matters, and even that cannot repair a proxy that does not measure p53 activity well."]),
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
<a href="#results">4. Results (4.1-4.8)</a> &middot;
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
TP53 wild-type cells, which provides a biological negative control on real data (Section 4.6).</p>
<div class="note"><b>Opening the loop removes one arrow, not the whole pathway.</b> Nutlin blocks the
reverse inhibitory arm MDM2 protein &#8866; p53, but it does not block p53's transcription-factor
activity. The forward branches therefore remain: p53 &#8594; MDM2 mRNA and p53 &#8594; CDKN1A mRNA.
Under Nutlin the observable target-gene structure is consequently the common-driver fork
MDM2 mRNA &#8592; p53 protein &#8594; CDKN1A mRNA, even though the negative-feedback cycle is open.</div>

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
<p>Unless a result states otherwise, the four-species simulations start from the basal state
(TP53 mRNA, p53 protein, MDM2 mRNA, MDM2 protein) = (15, 30, 20, 50). Figure 4 is the explicit
exception: to reproduce the start-up visualization in Figure 3 of Raharinirina et al. (2021), all
four species begin at zero and the mean-expression path therefore starts at the origin.</p>

<h4>Three-gene p53-MDM2-CDKN1A extension used in Figure 12D-E</h4>
<p>To test whether conditioning can separate a common-driver correlation from a direct regulatory
edge, we extend the four-species loop to six species: TP53 mRNA (m<sub>P</sub>), p53 protein (P),
MDM2 mRNA (m<sub>M</sub>), MDM2 protein (M), CDKN1A mRNA (m<sub>C</sub>) and CDKN1A protein (C).
The model has 13 reaction channels: transcription, mRNA degradation, translation and protein
degradation for each gene, plus MDM2-dependent p53 degradation. MDM2 and CDKN1A transcription
are separate Hill functions of the same p53-protein driver:</p>
<div class="formula">a<sub>MDM2 txn</sub>(P) = k<sub>M</sub> P^nM / (K<sub>M</sub>^nM + P^nM) + k<sub>M0</sub>
a<sub>CDKN1A txn</sub>(P) = k<sub>C</sub> P^nC / (K<sub>C</sub>^nC + P^nC) + k<sub>C0</sub>
a<sub>MDM2-dep p53 deg</sub>(P,M) = (1 - u) k<sub>MP</sub> M P / (K<sub>m</sub> + P)</div>
<p>Here u is Nutlin efficacy. Figure 12D-E uses u = 1, so the last propensity is zero: the
MDM2-protein &#8866; p53 arm is fully removed, while both p53-driven transcriptional arms remain.</p>
<table>
<tr><th>Process</th><th>Parameters used (per min; molecule-count constants where applicable)</th></tr>
<tr><td>TP53</td><td>transcription 3.0; mRNA decay 0.20; translation 4.0; basal p53-protein decay 0.02</td></tr>
<tr><td>MDM2 &#8866; p53</td><td>maximum MDM2-dependent p53 degradation 2.0; Michaelis constant K<sub>m</sub> = 120</td></tr>
<tr><td>p53 &#8594; MDM2</td><td>maximum transcription 6.0; K<sub>M</sub> = 80; Hill coefficient n<sub>M</sub> = 4; basal transcription 0.05</td></tr>
<tr><td>MDM2 turnover</td><td>mRNA decay 0.10; translation 1.2; protein decay 0.30</td></tr>
<tr><td>p53 &#8594; CDKN1A</td><td>maximum transcription 5.0; K<sub>C</sub> = 70; Hill coefficient n<sub>C</sub> = 3; basal transcription 0.05</td></tr>
<tr><td>CDKN1A turnover</td><td>mRNA decay 0.08; translation 1.0; protein decay 0.15</td></tr>
</table>
<p>We simulate <b>4000 independent cells</b> with the Gillespie direct algorithm, using initial state
(m<sub>P</sub>, P, m<sub>M</sub>, M, m<sub>C</sub>, C) = (15, 30, 20, 50, 15, 40) molecules,
random seed 93, and 21 retained observation times from 0 to 80 min in 4-min increments. These are
exact stochastic reaction counts with intrinsic molecular noise; no dropout, library-size factor,
negative-binomial sampling or other technical scRNA-seq noise is added. The full 4000 x 6 x 21
array is cached in <code>data/cache_step9_three_gene_nutlin.npz</code>, allowing Figures 12D-E to be
redrawn without rerunning the simulation. The model, accelerated simulator and cache-generation
code are saved in <code>src/model_three_gene.py</code>, <code>src/fast_ssa_three_gene.py</code> and
<code>analysis/step9b_three_gene_nutlin.py</code>, respectively. The cache also stores the species
order, initial state, parameter names and values, Nutlin efficacy, seeds and number of cells.</p>
<p>For Figure 12E, one time is sampled uniformly from the 20 post-Nutlin observations (4-80 min)
for each cell (snapshot seed 101). The two conditioning variables summarize the p53-protein input
seen by each target before measurement:</p>
<div class="formula">Z<sub>g</sub>(t) = sum_{{s&le;t}} H<sub>g</sub>(P(s)) exp[-gamma<sub>g</sub>(t-s)] Delta t,
H<sub>g</sub>(P) = P^ng / (K<sub>g</sub>^ng + P^ng),     g in {{MDM2, CDKN1A}}</div>
<p>The exponential weights account for loss of earlier transcriptional input through target-mRNA
decay. Partial correlation is then the Pearson correlation between the residuals after separately
regressing MDM2 mRNA and CDKN1A mRNA on the two-column matrix Z. This conditions on the actual
protein-level common driver available in the synthetic model. The noisy Figure 12E bars apply the
Splatter-style model of Section 3.5 to the two target transcript counts using seed 202, library-size
log-normal sigma 0.35, biological coefficient of variation 0.40, dropout midpoint 1.0 and dropout
shape -1.0. This produces 5.725% zero MDM2 observations, 5.700% zero CDKN1A observations, 10.875%
of cells with at least one zero target, and 0.55% of cells with both targets zero.</p>

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
<b>spurious</b> gene-gene correlation (Zappia et al., 2017); see Section 4.5.</div>

<h3>3.6 Direction and conditioning</h3>
<div class="formula">Granger (1-step): does A(t) improve prediction of B(t+1) beyond B(t)?
Partial correlation of B and C given A (dependence NOT explained by A):
  residual form:  rB = B - fit(B~A),  rC = C - fit(C~A);   pcorr(B,C|A) = corr(rB, rC)
  closed form:    pcorr(B,C|A) = ( rBC - rBA rCA ) / sqrt( (1 - rBA^2)(1 - rCA^2) )
     where rBC, rBA, rCA are ordinary Pearson correlations</div>
<p>Intuition: regress out the common driver A from both variables and correlate what remains; if B
and C are associated only because both depend on A, the residuals are uncorrelated (pcorr &asymp; 0).
Conditional dependence measures are nonlinear analogues (condition on A instead of removing a linear fit).
Granger causality is due to Granger (1969); distance
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
footing (Sections 4.5-4.6).</p>

{section_results()}

<h2 id="takeaways">5. Discussion and conclusions</h2>
<div class="key">
<ol>
<li>Intrinsic noise converts the deterministically-damped p53-MDM2 loop into sustained p53 pulses;
Nutlin-3 opens the loop, so p53 accumulates and oscillations cease (Sections 4.1-4.2).</li>
<li>The regulatory signal resides in higher-order moments (covariance, skewness), not in the mean of
the p53 transcript (Sections 4.3-4.4), as argued by Raharinirina et al. (2021).</li>
<li>Realistic technical noise (library size) and cell-cycle covariation create spurious
correlations that powerful dependence measures (Pearson, HSIC, distance correlation)
all flag as false positives; only the low-power kNN mutual-information estimator stays near zero.
The proper remedy is library-size normalization, not the choice of statistic (Sections 4.5-4.6).</li>
<li>On experimental MIX-seq data the p53 program (MDM2-CDKN1A co-regulation) is detected only in
TP53 wild-type cells under idasanutlin, with a clean isogenic mutant control (Section 4.6).</li>
<li>Symmetric measures report association only; the moment-based inference of Raharinirina et al.
(2021), and time-resolved/conditioned statistics, recover the directed edge p53&#8594;MDM2 and
distinguish regulation from mere correlation (Sections 4.7-4.8). Across five independent MBI
simulations, the p53&#8594;MDM2 direction remains positive in both closed-loop and Nutlin conditions,
whereas the exact edge weight varies; the robust conclusion is direction, not calibrated magnitude.</li>
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
