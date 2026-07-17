# Replicating minute-scale oscillatory sequences in entorhinal cortex

Replication of **Gonzalo Cogno et al. (2024), _Minute-scale oscillatory sequences in
medial entorhinal cortex_, Nature 625:338–344**, plus targeted extensions on public
data.

Prompted by [dandi/helpdesk#156](https://github.com/dandi/helpdesk/discussions/156).

**This is no longer only a replication.** The replication proper is the positive
control: the pipeline recovers the effect on the paper's own wheel-in-darkness
recordings. Everything past that is extension, and the space of possible extensions
is enormous — one could sweep every brain area and every condition on DANDI looking
for oscillatory sequences. That is not what this is. The extensions here are
*targeted*, each aimed at a question the paper itself raises and cannot answer with
its own data:

| extension | the paper's own words |
|---|---|
| navigation, free exploration, salient visual feedback (000053, 001701, 000690, EBRAINS open field) | "It remains an open question whether the ultraslow oscillatory sequences are present across a broader spectrum of behaviours, including sleep and free exploration, and in the presence of salient visual feedback" |
| task-free rest (000552) | same — the "sleep and free exploration" half |
| **is the oscillation just arousal?** (000690 pupil) | "consistent with a role for ascending neuromodulatory **arousal**-associated brain-stem circuits in controlling these oscillations" — an alternative the paper names but has no pupil, vascular or arousal measure to test |

## The claim being tested

MEC activity can organize into **ultraslow oscillations** (periods of tens of
seconds to minutes), and during those oscillations, cells fire in **periodic
sequences**. The effect was observed in a *sensory-minimized* condition: mice
head-fixed on a **rotating wheel in darkness**, 30–60 min, no scheduled rewards,
no change in location or heading.

**Important — darkness is where the effect was *observed*, not a condition shown
to be *required*.** Darkness/wheel was the *only* condition the paper ever ran;
there is no light or navigation comparison in it. The authors are explicit that
this is open (Discussion): *"It remains an open question whether the ultraslow
oscillatory sequences are present across a broader spectrum of behaviours,
including sleep and free exploration, and in the presence of salient visual
feedback. If so, it is possible that the sequences reset in the presence of
strong landmarks or sensory stimulation."*

So testing navigation datasets is not a defective replication — it is a partial
**answer to the paper's own stated open question**: do the sequences survive
navigation and salient visual input, and do they appear outside the mouse?

The paper's own specificity control was **region, not condition**: parasubiculum
and visual cortex (Fig. 5). And critically, it reports that *"while the calcium
activity of a fraction of cells in both brain areas was ultraslow and periodic,
in neither brain region were these oscillations organized into oscillatory
sequences."* Single-cell rhythmicity is non-specific **in the paper too** — only
the population sequence test discriminates. Our analysis reproduces that logic.

## Datasets

| Source | Species / region | Units | Duration | Role |
|---|---|---|---|---|
| **EBRAINS** [10.25493/SKKX-4W3](https://doi.org/10.25493/SKKX-4W3) | mouse MEC, Neuropixels 2.0 | 404 | 26 min | **positive control** — the original paper's own wheel-in-darkness data |
| ″ same file, `OpenField` epoch | mouse MEC, *same units* | 438 | 30 min | free foraging — **explicitly NOT analyzed in the paper**; an untested condition, not a validated negative |
| **EBRAINS two-photon**, mouse 60584 session 7 | mouse MEC, GCaMP6m | 522 | 60 min | **the paper's PRIMARY data** — its own Fig. 1b/2a–e example session |
| **DANDI 000053** (Mallory/Giocomo) | mouse MEC, Neuropixels | 74–408 | 26–44 min | VR straight track, rewards |
| **DANDI 001701** (Aery Jones 2026) | mouse MEC dorsal | ~150 | 40 min | X-maze navigation |
| **DANDI 000897** (Neupane/Fiete/Jazayeri) | **macaque** EC | 59 | 318 min | mental navigation — cross-species |
| **DANDI 000690** (Allen OpenScope) | mouse **MEC + PaS + V1**, simultaneous | up to 726 MEC | 120–135 min | **the paper's Fig. 5 design within a session**; passive viewing, no darkness |
| **DANDI 000552** (Huszár et al.) | mouse **CA1** | 102–221 | 39–148 min | **task-free rest** — region control in the one condition the archive can supply |

### Screened and rejected

Full catalogue in **[SCREENING.md](SCREENING.md)** — every entorhinal or long-rest
candidate on DANDI, from a scan of all 875 dandisets, with the quantitative reason
each failed. Two structural findings drive most of it:

- **The 17.6-minute floor.** The paper's Welch window is 17.6 min (8,192 × 129 ms),
  so one analysis window needs that much *continuous* recording in one condition.
  It is what disqualifies IBL (its rest block is capped at exactly 10.0 min by
  protocol), 000053 (8.4 min of darkness) and 000059 (12.25 min of septal cooling).
- **The condition is never in the metadata.** Of 875 dandisets, 3 mention darkness
  (none relevant) and zero mention ultraslow/infraslow oscillation. Rest and dark
  blocks are found only by streaming files — and the derived data usually stops at
  the task boundary anyway (001634 has ten ideal 17–33 min rest epochs and sorted
  units in *none* of them).

The EBRAINS **wheel** epoch is the methodological keystone: a detector that
cannot fire there cannot be trusted to report a null anywhere else.

The `OpenField` epoch in the same file is *not* a validated negative control. The
paper's Methods state these trials "were not used in the present study", so the
effect there is simply **unknown**. Our result on it (no population sequences,
same units, same probe) is therefore a *new observation* bearing on the authors'
open question about free exploration — not a check on the pipeline.

## Method (per Methods of the paper, ephys variant)

1. bin spikes at **120 ms**
2. convolve with a **Gaussian kernel, σ = 5 s**
3. **binarize** at mean + 1 SD
4. autocorrelation + **power spectral density**; an ultraslow oscillation is a
   PSD peak **below 0.1 Hz** absent from shuffled data
5. **PCA** on the (time × cell) activity matrix; each cell's angle **θ** = angle
   of its (PC1, PC2) loading vector; sorting by θ exposes periodic sequences

### Two statistical traps found and fixed

Both were caught *because* the positive control existed — each produced a clean,
confident, and completely wrong null result on data where the effect is known to
be present.

**1. Band power cannot detect the effect.** The σ = 5 s kernel is a low-pass
filter, so ~75% of the power lands in the ultraslow band *for the real data and
for the shuffled surrogates alike* (measured: 0.74 vs 0.77). A statistic based
on power-in-band is blind by construction. The signature is a **narrowband peak
above the broadband 1/f-like background**, so the statistic must whiten that
background: `max over band of PSD_data(f) / mean_surrogate_PSD(f)`.

**2. The permutation resolution floor makes FDR unreachable.** With 200
surrogates the smallest attainable p-value is 1/201 ≈ 0.005, but
Benjamini–Hochberg across 404 cells demands p < 0.05/404 ≈ 1.2e-4. *No cell can
ever be significant*, no matter how strong the oscillation — the original data
returned 0/404 with a max excess-power ratio of 47×. Fixed by pooling the null
across cells (the normalized statistic is exchangeable under H0), giving ~1e-5
resolution.

The choice of which signal to run the spectrum on (binarized activity vs.
smoothed rate vs. spike-train autocorrelogram) is settled empirically in
`03_detector_bakeoff.py` against the **wheel** epoch — the one condition with a
known-present effect — not by intuition. (That script's stated goal of "low on
open field" was mistaken; the open field is an untested condition, not a known
negative. The detector was selected on the wheel column regardless: only the
autocorrelogram spectrum recovers the effect, 64% vs 0–2%.)

### The population sequence test (the load-bearing criterion)

Single-cell rhythmicity is non-specific, so the actual replication test is at the
population level, with two design choices:

- **Null = independent circular shift of each cell's trace.** A circular shift
  leaves every cell's own spectrum and autocorrelation exactly intact and destroys
  only the timing *between* cells — so it isolates the paper's real claim (a shared
  oscillation with cells at staggered phases) from "cells are individually slow".
- **Windowed, not whole-session.** A whole-session statistic averages the sequence
  epochs together with quiet stretches and washes out intermittent sequences — it
  gave a false negative on the paper's own mouse 102335. We slide a 300 s window,
  score each against the shift null, and take a binomial test on the count of
  significant windows (the "oscillation score"). This is the paper's own
  fraction-of-session-oscillating logic, with a per-window null attached, and it
  fires on both original Neuropixels mice.

> **Note.** The windowed binomial as originally parameterized (300 s window, 100 s
> step) counts **overlapping** windows and is anti-conservative. All numbers in the
> results and interpretation below are the original overlapping values, marked
> **[†](#the-overlapping-window-caveat)**. They inflate the positives — the 2/2
> positive control is really 1/2 — but change **no conclusion**, because the nulls
> are unaffected. See **[The overlapping-window caveat](#the-overlapping-window-caveat)**
> for why, and for the corrected numbers side by side.

## Results (every session, every dataset)

The paper's reported quantity is the **fraction of sessions with periodic
sequences** (Fig. 5g: MEC 15/27, PaS 0/25, VIS 0/19), because the effect is
strongly session-variable. We compute the same with the windowed sequence test
(`figures/SUMMARY_all_sessions.png`). Numbers are the original overlapping-window
values, marked **[†](#the-overlapping-window-caveat)**; the corrected values are in
[that section](#the-overlapping-window-caveat).

| group | condition | sessions with sequences **[†](#the-overlapping-window-caveat)** |
|---|---|---|
| **EBRAINS wheel / darkness** | positive control | **2/2 = 100%** |
| EBRAINS open field (same units) | untested in paper | 0/2 |
| **DANDI 000053** mouse MEC | VR straight track | **0/20 = 0%** |
| **DANDI 001701** mouse MEC | X-maze | **3/114 = 3%** |
| DANDI 001701 visual cortex | region control | 9/110 = 8% *(confound)* |
| DANDI 000897 macaque EC | mental navigation | 8/15 = 53% *(confound)* |
| **DANDI 000690** mouse MEC | passive viewing | 3/11 = 27% |
| DANDI 000690 visual cortex | *simultaneous* control | 1/11 = 9% |
| DANDI 000690 parasubiculum | *simultaneous* control | 1/3 |
| **DANDI 000552** mouse CA1 | task-free rest | **0/6 = 0%** |

(The **two-photon** primary-data result is reported separately in
"[The paper's primary data](#the-papers-primary-data)".)

**The pipeline reproduces the original effect.** On the paper's own wheel-in-darkness
data it recovers ultraslow rhythmicity in 64% of MEC cells (paper: ~91% by a
different criterion), a median period of 39 s, and coherent population sequences in
**both** Neuropixels mice (104638 p≈1e-6; 102335 p=0.025), while the two open-field
epochs from the same units are negative. (Under the corrected test the positive
control is 1/2 — 104638 survives, 102335 does not; see
[the caveat](#the-overlapping-window-caveat).)

**We do not replicate minute-scale oscillatory *sequences* in mouse MEC during
navigation.** Across **134 mouse MEC navigation sessions** (0/20 VR track + 3/114
X-maze) the sequence rate is at chance — and the within-session visual-cortex
region control is *higher* (9/110), not lower, which is the opposite of an
intrinsic MEC-specific effect.

### Two results that matter more than the bare null

**Single-cell ultraslow rhythmicity is not specific and must not be the
replication criterion.** It is high *everywhere*: ~50–75% of cells in every
group, including open field and visual cortex. A grid/place cell sweeping its
fields is slowly rate-modulated for trivial spatial reasons — and the paper's own
Fig. 5 says the same (PaS and VIS cells are "ultraslow and periodic" yet form no
sequences). Any replication that stops at single-cell spectra reports a false
positive. The load-bearing test is the *population* one.

**Behaviour/task structure manufactures the population signature — in two species.**
The two apparent "positives" outside the wheel are both behavioural artifacts,
each nailed by a dedicated control:

- *Mouse visual cortex* (9/110, higher than MEC). The animal runs laps, the
  visual scene repeats, V1 tiles the lap cycle. `05_behavior_control.py`
  (`figures/behavior_control_001701.png`): position, MEC-PC1 and V1-PC1 all peak
  at 0.007–0.010 Hz — the same lap rhythm.
- *Macaque EC* (8/15, at first glance a striking cross-species replication).
  `10_macaque_control.py` (`figures/macaque_task_control.png`): the population PC1
  tracks trial-onset density (it drops to zero during task pauses and resumes with
  the trials), and the PC1 spectrum overlaps the trial-onset spectrum across the
  whole band (|corr| 0.30–0.53). It is task-engagement block structure over a
  multi-hour session, not an intrinsic rhythm.

### CA1 at rest: single-cell rhythmicity is meaningless

`11_ca1_rest_000552.py`, `figures/ca1_rest_000552.png`. **0/6 sessions with
sequences** (102–221 units, 39–148 min of continuous task-free awake rest; median
oscillation score 0.00, Wilcoxon vs 0.05 chance p=0.031).

The null itself is *expected* — the paper's contrast is that MEC carries a
low-dimensional code with ~94% of cells locked to the sequence whereas hippocampal
sequences involve ~5% of the network. The **dissociation** is the point, and it is
more extreme here than anywhere else in the panel:

| | single-cell rhythmicity | sessions with sequences |
|---|---|---|
| CA1, task-free rest | **83%** (highest in the panel) | **0/6** |
| ORIGINAL wheel / darkness | 53% | 2/2 **[†](#the-overlapping-window-caveat)** |

CA1 at rest is *more* ultraslow-rhythmic than the paper's own data where the effect
is known present, at a median peak of 0.0066 Hz — inside the paper's reported
0.006–0.008 Hz band — with no population sequences at all. At its best 500 s window
it reaches rotation 0.22, against 0.75 for the wheel control. A pipeline that stopped
at single-cell spectra would rank CA1-at-rest as a *stronger* replication than the
original data.

#### Not a power artifact — tested directly

`13_power_analysis.py`, `figures/power_analysis.png`. CA1 rest has 1.9× fewer units
than the wheel control (215 vs 404) and 2.5× lower firing rates (1.25 vs 3.18 Hz), so
"too few / too sparse cells" was a live alternative. We degraded the positive control
to CA1's conditions and re-ran the identical test:

| positive control, degraded | sequences detected |
|---|---|
| n = 50 units (a quarter of CA1's) | **5/5** |
| **n = 215 AND thinned to 1.25 Hz (fully CA1-matched)** | **5/5 in both mice** |
| *CA1 rest, actual* | *0/6* |

The effect survives every handicap CA1 imposes, and survives even at n=50 —
unsurprising given ~94% of MEC cells participate; a population-wide effect does not
need a large sample. Insufficient power is not an available explanation. Region,
species and state remain uncontrolled.

### MEC vs its own control regions, recorded simultaneously (000690)

`14_openscope_000690.py`, unit-matched within session. **MEC 3/11 (27%), V1 1/11
(9%), PaS 1/3 [†](#the-overlapping-window-caveat)**.

This is the only dataset that runs the paper's Fig. 5 contrast *within* a session, so
the stimulus is common to all regions by construction and cannot manufacture a
MEC-vs-V1 difference. But **the contrast does not reach significance** (Fisher exact
MEC vs V1, p=0.59 — true on either window scheme), and it does not survive inspection
of the individual sessions:

- **714615** is the one clean MEC-specific positive (survives the window correction:
  9/26 independent windows, p<1e-4, rotation 0.428 vs null 0.299).
- **702134** is positive in **MEC, V1 and PaS at once** — a global signal, not an MEC
  one. It is also the session with the strongest MEC–pupil coherence (z=+3.9), so
  this is most likely the arousal signal appearing in every region simultaneously.
- **695764** was a subsampling artifact: significant at the 250-unit subsample, gone
  (p=0.69) with all 401 units.

So **one** MEC session shows a defensible region-specific sequence, versus one for
V1 — the within-session contrast this dataset was added to provide does not hold.
MEC does clear the 5% chance level across sessions (binomial p=0.015) where V1 does
not (p=0.43), which is real but much weaker than "MEC beats its simultaneous control".

The rest of the panel still orders wheel/darkness → passive viewing → navigation,
consistent with the authors' speculation that "sequences reset in the presence of
strong landmarks or sensory stimulation". It is suggestive, not established: the
groups differ in lab, rig, species-strain and sorting as well as condition.

### Is the oscillation just arousal? (000690 pupil) — the paper's orphaned alternative

`15_arousal_control.py`, `figures/arousal_control_000690.png`. **Yes, substantially.**

The paper names ascending **arousal** circuits as the likely driver of the
single-cell oscillation and cites infraslow pupil oscillation and vascular confounds
— but has no pupil, vascular or arousal measure anywhere, and 0.006–0.008 Hz sits in
the arousal/vasomotion range. 000690 carries EyeTracking alongside the archive's best
MEC coverage, so the control is runnable.

Coupling is tested as band coherence between the MEC population and pupil area,
against an independent **circular shift** of the pupil trace — the same null as the
sequence test, preserving each signal's own spectrum and destroying only their
relative timing.

| | median across 11 sessions |
|---|---|
| MEC–pupil band coherence, observed | **0.187** |
| same, circular-shift null | 0.086 |
| coherence z vs null | **+3.9** |
| sessions significant at p<0.05 | **8/11** |
| \|peak frequency difference\| | 0.0033 Hz |

MEC population peaks land at 0.005–0.012 Hz — squarely in the paper's own
0.006–0.008 Hz band — and are coherent with pupil at ~2.2× the shift null in 8 of 11
sessions. **The ultraslow oscillation in MEC substantially tracks arousal**, which is
what the authors themselves proposed and could not test.

This does **not** touch the sequence claim. It supports the paper's own two-part
model: the oscillation is ascending and global (and now measurably arousal-coupled),
while sequence formation may still be MEC-specific. It does mean that any single-cell
ultraslow result — in MEC or anywhere else — should be read as an arousal measurement
until a population sequence test says otherwise.

A trap worth recording: correlating the two **log spectra** gives r≈0.94 and is
meaningless, because both signals are 1/f-like and any two such spectra agree in log
space whether or not they share an oscillation. That is the same failure as "band
power cannot detect the effect" above. Only the shift-null comparison is informative.

### Interpretation

**The replication holds, on both modalities.** Minute-scale oscillatory sequences
reproduce on the paper's own wheel-in-darkness Neuropixels data and — the part that
actually carries the paper's claims — on its **primary two-photon data**: session
60584/7 gives an oscillation score of 0.21 (p=0.0013 [†](#the-overlapping-window-caveat))
and a near-circular PC1–PC2 ring (ratio 1.15). Both survive the window correction
(104638 and the imaging session pass; only the paper's weaker mouse 102335 does not,
dropping the positive control to 1/2 — see [the caveat](#the-overlapping-window-caveat)).
The single-cell half of the imaging session reproduces independently at the reported
~0.0066 Hz
([rly/replicate-gonzalo-cogno-2023](https://github.com/rly/replicate-gonzalo-cogno-2023)).
A detector was validated on known-positive data before any null was trusted anywhere.

**Sequences are not detectable in mouse MEC during navigation** (0/20 VR track +
3/114 X-maze [†](#the-overlapping-window-caveat)), and only weakly under passive
viewing (3/11), consistent with the authors' speculation that "sequences reset in the
presence of strong landmarks or sensory stimulation". Across the panel the rate orders
wheel/darkness → passive viewing → navigation, which is suggestive but confounded by
lab, rig and sorting. (These nulls are unchanged by the window correction.)

**Single-cell ultraslow rhythmicity is not evidence for the claim, and is now
partly explained.** It is high everywhere (~50–90% in every group, highest of all in
CA1 at rest, which has zero sequences), and in MEC it is **coherent with pupil-indexed
arousal in 8/11 sessions** at ~2.2× a circular-shift null. That is the paper's own
proposed mechanism for the oscillation, tested for the first time. Any single-cell
ultraslow result should be read as an arousal measurement until a population sequence
test says otherwise.

**What survives as the paper's distinctive claim** is therefore the narrower and more
interesting one: not that MEC cells oscillate slowly — that appears to be arousal,
and it is everywhere — but that MEC alone organizes those oscillations into periodic
population *sequences*. Nothing here contradicts that; it reproduces on the primary
imaging data and on one Neuropixels control. But the within-session MEC-vs-V1 contrast
in 000690 does **not** confirm region specificity (3/11 vs 1/11, not significant; the
one all-region-positive session tracks pupil), so MEC-specificity rests on the
paper's own data plus the region *nulls* here, not on a positive out-of-sample
discrimination.

**What this does not license**: that darkness is *necessary*. The paper never ran
another condition, so necessity is unestablished, and a null cannot separate "sensory
drive abolishes the sequences" from "insufficient power" — though for CA1, power is
now ruled out directly.

Caveats: the detector is calibrated to the stronger wheel session (104638); weaker
intermittent sequences of the kind the paper notes in mouse 102335 could be missed.
The macaque has few units per session (25–64) and a trial-structured task. 000690
session 702134 is positive in all three regions at once, which is a global signal, not
an MEC-specific one.

Depositing the EBRAINS wheel data on DANDI would be valuable — it is the only public
recording in the condition where the effect is known to exist.

## The overlapping-window caveat

Every windowed `sequences: yes/no` call in this report (the **[†](#the-overlapping-window-caveat)**
marks) comes from `binomtest(n_significant_windows, n_windows, 0.05)` over a 300 s
window sliding in **100 s steps**. Adjacent windows therefore share **two-thirds of
their samples**, while the binomial assumes they are **independent**. So it counts
~3× more evidence than exists and is **anti-conservative** — it manufactures
significance.

**This is the repo's parameterization, not the paper's.** The paper does use a
binomial (Extended Data Fig. 8, the travelling-wave analysis), but only over
genuinely **independent** units — 15 distinct sessions, or non-overlapping spatial
bins (e.g. "1/15, probability = 0.37, binomial distribution"; "1/8, probability =
0.28"). Independence is the assumption every binomial in the paper satisfies and the
overlapping-window version violates. And the paper's own session-level call is not a
windowed binomial at all: its **oscillation score** is a single whole-session Welch
PSD of the population phase with a peak-prominence test (peak > 9× the tail mean and
> 9× the pre-peak minimum), and its **sequence score** is thresholded at the 99th
percentile of a per-session shuffle. The sliding-window binomial is this repo's
approximation of that, and the overlap is where it goes wrong.

`windowed_sequence_test` now also returns the statistic on **non-overlapping**
windows (every 3rd window of the same grid, so free to compute), and
`17_window_independence.py` re-runs the whole panel both ways
(`figures/window_independence.png`). The correction is **asymmetric**, which is why
no conclusion moves:

| group | overlapping (reported above) | **independent (corrected)** |
|---|---|---|
| EBRAINS wheel (positive control) | 2/2 = 100% | **1/2 = 50%** |
| EBRAINS open field | 0/2 | 0/2 |
| two-photon primary data | p=0.0013 (7/34) | **p=0.0022 (4/12), still YES** |
| 000053 MEC (VR track) | 0/20 | **0/20** |
| 001701 MEC (X-maze) | 3/114 = 3% | **2/114 = 2%** |
| 001701 V1 | 9/110 = 8% | **6/110 = 5%** |
| 000897 macaque EC | 8/15 = 53% | **6/15 = 40%** |
| 000690 MEC (passive) | 3/11 = 27% | **2/11 = 18%** |
| 000690 V1 | 1/11 = 9% | **1/11 = 9%** |
| 000552 CA1 (rest) | 0/6 | **0/6** |

- **Every null is unchanged** — a liberal test that found nothing finds nothing when
  made stricter. The CA1 and navigation stories, which are carried by nulls, do not
  move.
- **Every positive weakens** in proportion to how much it leaned on overlap.
- **The positive control drops to 1/2.** Mouse 102335 goes from 3/13 windows
  (p=0.025) to 1/5 (p=0.23) — and its median observed window rotation (0.498) is
  *below* its median null (0.588), so that "positive" was pure overlap. Mouse 104638
  survives cleanly (3/5, p=0.0012), as does the two-photon primary session. So the
  detector is still validated on known-positive data, on the *stronger* of the paper's
  two Neuropixels mice plus its own imaging session — just not on 102335, the subtle,
  subset-dependent one the paper itself flags.

Note that `09_validate_windowed.py` still states the detector "fires on both original
wheel mice"; under the independent test that is 1/2. The overlapping numbers are kept
as the headline only to report the replication *as the original pipeline computes it*;
the independent column is the statistically defensible one.

**Is 102335 "really" oscillatory by the paper's own method?** Not resolvable on public
code. The repo's windowed-rotation test is not the paper's classifier — the paper uses
an *oscillation score* (a population-phase PSD gate plus a τ–d joint-distribution score
over concatenated sequences). Both top-level functions are in the authors' repo, but
the pipeline cannot be run end-to-end: leaf dependencies are missing, including two
that set the detection thresholds, and the released scorer diverges from the Methods
text (`nperseg = 2048` not 8,192; a **4×** prominence threshold not 9×, commented
`%Used to be 9`). So whether 102335 is oscillatory is not settled by anything runnable
here; the repo finding no sequence in it is a **detector difference from the paper's
method, not evidence against the paper** (the one fully-specified component, the
population-phase PSD peak, is if anything *cleaner* for 102335 than for 104638). A
**future direction** is to port the paper's population-phase gate — fully specified in
the released code — and run it across the panel as a second, paper-native classifier
alongside the binomial.

## Where else to look

Candidate regions connected to MEC, and what data exists for each, are tabulated in
[SCREENING.md](SCREENING.md#connected-regions-candidates-and-their-data). The short
version: **lateral entorhinal cortex** is the strongest untested candidate — adjacent,
reciprocally connected, invoked by the paper itself for its slow drift, and never
recorded in the paper — and the only LEC dandiset is empty. **Medial septum** is the
way to ask whether the rhythm is imposed rather than intrinsic; all three septum
datasets fail as deposited. Neither gap is fixable by analysis.

### The paper's primary data

The paper's headline numbers — 91% of cells oscillatory, ~94% phase-locked, 44% at
0.006–0.008 Hz, the ring manifold, 15/27 oscillatory sessions — are **not** from the
Neuropixels recordings used as the positive control here. They are from **two-photon
imaging of 6,231 cells across 15 sessions**. That is the actual claim.

**The single-cell half was already replicated**, independently, in
[rly/replicate-gonzalo-cogno-2023](https://github.com/rly/replicate-gonzalo-cogno-2023):
it goes to the EBRAINS two-photon deposit (`data/calcium/60584/2019-01-29/MUnit_0`,
the paper's own Fig. 1 example session), applies the paper's preprocessing, and
reproduces **Figure 1** — the same vertical banding in the stacked autocorrelations,
example cells peaking at ~0.0066 Hz with harmonics, matching the reported value.

**The sequence half is tested here** (`16_imaging_sequences.py`), which is the part
that matters: single-cell rhythmicity is exactly what this repo shows is non-specific.
Result on the same session, 522 cells, 60 min, using the paper's own activity matrix:

| | |
|---|---|
| significant windows | **7/34 [†](#the-overlapping-window-caveat)** |
| oscillation score | **0.21** |
| p_session | **0.0013** |
| **sequences** | **YES** |
| ring geometry (PC1/PC2 ratio) | **1.15** (a ring is ~1.0) |

For reference the wheel Neuropixels sessions score 0.54 (p=1e-6) and 0.23 (p=0.025),
so the imaging session sits alongside the weaker mouse. And unlike that weaker mouse
(102335, which fails the [window correction](#the-overlapping-window-caveat)), **the
imaging session survives it** — 4/12 non-overlapping windows, p=0.0022 — so this
result does not depend on the flaw. The PCA-sorted raster reproduces Fig. 2b and the
PC1–PC2 manifold reproduces the Fig. 2c ring at a near-circular 1.15 aspect
(`figures/imaging_sequences.png`).

So the paper's distinctive claim now stands on its primary data, not only on two
Neuropixels mice.

#### A third statistical trap, and this one nearly produced a false negative

Applying the ephys pipeline's convention to the imaging data — Gaussian σ = 5 s, then
z-score — returns **0/34 significant windows, p=1, "no sequences"** on data where the
effect is present. The reason is instructive: deconvolved calcium is already
temporally smooth, so smoothing again oversmooths, and the PC1/PC2 trajectory then
advances consistently *for the circular-shift surrogates too*. Observed rotation is a
high-looking **0.853 — against a null of ~0.83 (z=+1.2)**. The statistic saturates and
loses all discriminative power. The σ = 5 s kernel exists in the ephys path to mimic
calcium dynamics that imaging data already has; applying it to calcium double-counts
it, which is why the paper binarizes imaging data instead.

Smoothing is right for the manifold *picture* (the paper smooths for Fig. 2c) and
wrong for the *statistic*. Same family as the two traps above: a confident, clean,
completely wrong answer, caught only because the effect was known to be present.

#### Two reproducibility findings from the imaging work

- **The authors' MATLAB cannot be run as published.** Both `analysis_npx_b.m` and
  `Autocorrelations_PSDs_examples_b.m` call `npx_init`, which is not in the
  repository, and paths are hardcoded to a lab drive (`W:\npxwaves\MEC\`). So
  *computational* reproducibility — same code, same data, same numbers — is still
  blocked, independent of the methodological reproduction here.
- **Cell-count discrepancies.** The paper reports 469 good units for mouse 104638; the
  deposited file yields 487. It reports 484 cells for imaging session 60584/7; `iscell`
  gives 522, and the stated SNR>4 criterion is non-binding (reimplemented on
  F − 0.7·Fneu it gives a median SNR of 21.5, with no cell falling in 0 < SNR ≤ 4).
  Neither gap is material — ~94% of cells participate — but neither reproduces exactly.

### A methods bug worth flagging for anyone reusing this

001701 labels brain regions under **two naming conventions across subjects** —
spelled-out ("Entorhinal area medial part dorsal zone") for some, Allen atlas
codes ("ENTm2/ENTm3/ENTm5", "VISpl5") for others — and the units→electrode map
must be read from each unit's electrode **sub-table**, not by indexing a location
array with electrode *ids* (which are non-contiguous on multi-shank probes). Doing
either wrong silently drops or mis-assigns whole subjects. The corrected loader
finds MEC units in 114/218 sessions across 7 subjects; the other ~104 sessions are
hippocampal recordings with no sorted MEC units. See `scan_001701_regions.py`.

## Running it

```bash
pip install -r requirements.txt

# resumable run over every session x region (writes results/sessions/*.json)
python3 scripts/07_run_all_sessions.py

# aggregate to the fraction-of-sessions figure the paper reports
python3 scripts/08_session_summary.py

# CA1 during long task-free rest (000552) + its dissociation figure
python3 scripts/11_ca1_rest_000552.py
python3 scripts/12_ca1_rest_figure.py

# is the CA1 null just insufficient power? (subsample + rate-match the control)
python3 scripts/13_power_analysis.py

# MEC vs its own control regions, simultaneous (000690); then the arousal control
python3 scripts/14_openscope_000690.py
python3 scripts/15_arousal_control.py

# the sequence test on the paper's primary two-photon data (downloads from EBRAINS)
python3 scripts/16_imaging_sequences.py
```

All data is streamed from DANDI and EBRAINS on demand (remfile + a local disk
cache under `cache/`, which is git-ignored); nothing is downloaded in bulk. The
run is resumable, so it can be stopped and restarted freely. Paths resolve
relative to the repository, so no configuration is needed.

## Layout

```
scripts/
  00_triage_dandisets.py   # screen DANDI for entorhinal datasets
  01_inspect_sessions.py   # stream units tables; units, duration, regions
  mec_oscillation.py       # core: activity matrices, PSD, surrogates, PCA sequences
  loaders.py               # per-dataset NWB quirks + EBRAINS control
  diag_psd.py              # shows why band power cannot detect the effect (trap 1)
  03_detector_bakeoff.py   # pick the detector using wheel-vs-openfield ground truth
  04_analyze.py            # the validated pipeline over the whole panel -> figures
  05_behavior_control.py   # mouse V1: is the 'oscillation' just the lap cycle?
  06_summary_figure.py     # single-session cross-dataset summary
  07_run_all_sessions.py   # RESUMABLE run over every session x region -> results/sessions/
  08_session_summary.py    # fraction of sessions with sequences (paper's quantity)
  09_validate_windowed.py  # windowed detector validated on the original wheel mice
  17_window_independence.py# re-run the panel with non-overlapping windows (the correction)
  10_macaque_control.py    # macaque EC: is the 'oscillation' task-engagement structure?
  11_ca1_rest_000552.py    # CA1 during long TASK-FREE rest (the condition, minus the region)
  12_ca1_rest_figure.py    # the rhythmicity-vs-sequences dissociation + rasters
  13_power_analysis.py     # degrade the positive control to CA1's n and rate: does it still fire?
  14_openscope_000690.py   # MEC vs PaS vs V1 recorded simultaneously (unit-matched)
  15_arousal_control.py    # is the ultraslow oscillation just pupil-indexed arousal?
  16_imaging_sequences.py  # the sequence test on the paper's PRIMARY two-photon data
  scan_001701_regions.py   # metadata scan: which sessions actually have MEC/V1 units
figures/  results/  cache/
SCREENING.md                       # every dandiset screened, and why it was rejected
replication_notebook.py / .ipynb   # narrative walkthrough
```
