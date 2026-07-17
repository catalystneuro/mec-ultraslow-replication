# Replicating minute-scale oscillatory sequences in entorhinal cortex

Replication attempt of **Gonzalo Cogno et al. (2024), _Minute-scale oscillatory
sequences in medial entorhinal cortex_, Nature 625:338–344**, extended to public
datasets other than the original.

Prompted by [dandi/helpdesk#156](https://github.com/dandi/helpdesk/discussions/156).

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
| **DANDI 000053** (Mallory/Giocomo) | mouse MEC, Neuropixels | 74–408 | 26–44 min | VR straight track, rewards |
| **DANDI 001701** (Aery Jones 2026) | mouse MEC dorsal | ~150 | 40 min | X-maze navigation |
| **DANDI 000897** (Neupane/Fiete/Jazayeri) | **macaque** EC | 59 | 318 min | mental navigation — cross-species |
| **DANDI 000552** (Huszár et al.) | mouse **CA1** | 102–221 | 39–148 min | **task-free rest** — region control in the one condition the archive can supply |

Screened and rejected: **000582** (only 3–5 units/session), **000638** (190-min
sessions but too few MEC units), **000943** (raw ecephys only — no spike sorting).

**000552 buys the condition by giving up the region.** Every entorhinal dataset on
DANDI runs a task from start to finish, so none can be asked the paper's question
in the paper's condition (no task, no rewards, minimal sensory drive). 000552 is
the archive's best approximation: multi-hour sessions where unlabelled PRE/POST
home-cage blocks bracket a maze epoch. 24 of 64 sessions have ≥100 units and ≥30
min of *continuous* task-free awake rest. The epochs table carries only
start/stop times, so task-free blocks are identified empirically as those
containing no trials, and brain state comes from `processing/behavior/SleepStates`.
Bouts are never concatenated — that would destroy the phase continuity a
minute-scale sequence lives in. Caveats: `electrodes.location` is `unknown`
throughout (CA1 comes from the publication), and darkness is undocumented.

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
  fires on **both** original Neuropixels mice. The detector is calibrated to the
  strong wheel session and deliberately not tuned to force borderline cases.

## Results (every session, every dataset)

The paper's reported quantity is the **fraction of sessions with periodic
sequences** (Fig. 5g: MEC 15/27, PaS 0/25, VIS 0/19), because the effect is
strongly session-variable. We compute the same, using the windowed sequence test
(`figures/SUMMARY_all_sessions.png`):

| group | condition | sessions with sequences |
|---|---|---|
| **EBRAINS wheel / darkness** | positive control | **2/2 = 100%** |
| EBRAINS open field (same units) | untested in paper | 0/2 |
| **DANDI 000053** mouse MEC | VR straight track | **0/20 = 0%** |
| **DANDI 001701** mouse MEC | X-maze | **3/114 = 3%** |
| DANDI 001701 visual cortex | region control | 9/110 = 8% *(confound)* |
| DANDI 000897 macaque EC | mental navigation | 8/15 = 53% *(confound)* |
| **DANDI 000552** mouse CA1 | **task-free rest** | **0/6 = 0%** |

**The pipeline reproduces the original effect.** On the paper's own
wheel-in-darkness data it recovers ultraslow rhythmicity in 64% of MEC cells
(paper: ~91% by a different criterion), a median period of 39 s, and coherent
population sequences in **both** Neuropixels mice (104638 p≈1e-6; 102335 p=0.025),
while the two open-field epochs from the same units are negative.

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

### CA1 at rest: the sharpest evidence that single-cell rhythmicity is meaningless

`11_ca1_rest_000552.py`, `figures/ca1_rest_000552.png`. **0/6 sessions have
sequences** (102–221 units, 39–148 min of continuous task-free awake rest).

The value of this group is not the null itself — a CA1 null is *expected*, since
the paper's contrast is that MEC carries a low-dimensional code with ~94% of cells
locked to the sequence whereas hippocampal sequences involve ~5% of the network.
It is the **dissociation**, which is more extreme here than anywhere else in the
panel:

| | single-cell rhythmicity | sessions with sequences |
|---|---|---|
| CA1, task-free rest | **83%** (highest in the panel) | **0/6** |
| ORIGINAL wheel / darkness | 53% | 2/2 |

**CA1 at rest is *more* ultraslow-rhythmic than the paper's own data where the
effect is known to be present, at a median peak of 0.0066 Hz — squarely inside the
paper's reported 0.006–0.008 Hz band — and yet has no population sequences at
all.** Median oscillation score is 0.00 (Wilcoxon vs the 0.05 chance level,
p=0.031). At its single best 500 s window CA1 reaches a rotation index of 0.22,
against 0.75 for the wheel control.

This is what the paper's own Fig. 5 asserts (PaS and VIS cells are "ultraslow and
periodic" yet form no sequences), reproduced in a third region and in a condition
the authors never tested. It also sets a trap for anyone reusing this work: a
pipeline that stops at single-cell spectra would report CA1-at-rest as a *stronger*
replication than the original data. Only the population test discriminates.

#### Is the CA1 null just too few / too sparse cells? No — tested directly.

`13_power_analysis.py`, `figures/power_analysis.png`. The worry is real on its
face: CA1 rest has **1.9× fewer units** than the wheel control (215 vs 404) and
**2.5× lower firing rates** (1.25 vs 3.18 Hz median). Either could hide a real
effect. So we degraded the positive control — the paper's own wheel data, effect
known present — down to CA1's conditions and re-ran the identical test.

| positive control, degraded | sequences detected |
|---|---|
| n = 50 units (a quarter of CA1's) | **5/5** |
| n = 215 units (CA1's count) | 4/5 and 5/5 |
| **n = 215 AND thinned to 1.25 Hz (fully CA1-matched)** | **5/5 in both mice** |
| *CA1 rest, actual* | *0/6* |

**The effect survives every handicap CA1 imposes.** Subsampled to CA1's exact unit
count and thinned to CA1's exact firing rate, both original wheel sessions still
show sequences in 5/5 draws (median oscillation score 0.38 and 0.31, against CA1's
0.00). It survives even at n=50 — unsurprising given the paper reports ~94% of MEC
cells locking to the sequence: a population-wide effect does not need a large
sample. Thinning is the right SNR control because thinning a rate-modulated
process preserves the modulation and adds Poisson noise.

So insufficient power is **not** an available explanation for the CA1 null. The
null stands on its own.

What remains uncontrolled is everything else that differs: region (CA1 vs MEC),
species (rat vs mouse), and state (home-cage rest vs running a wheel in darkness).
This is not the paper's condition, and a null still cannot prove MEC-specificity —
but "not enough cells" is now ruled out rather than merely doubted.

### Interpretation

Minute-scale oscillatory **sequences** reproduce on the original wheel data and
are **not detectable in mouse MEC during navigation** (0/20 + 3/114 sessions),
consistent with the authors' own speculation that "sequences reset in the presence
of strong landmarks or sensory stimulation". This is a partial answer to the open
question they pose. The apparent positives in mouse V1 and macaque EC are
behavioural/task confounds, not evidence for the phenomenon.

What this does **not** license is the claim that darkness is *necessary* — the
paper never tested another condition, so necessity is unestablished, and a null
cannot separate "sensory drive abolishes the sequences" from "insufficient power".
But the multi-session scale substantially strengthens the null: 134 MEC navigation
sessions, with a validated detector that fires on both original wheel sessions and
is deliberately *not* tuned to force borderline cases.

Caveats: the detector is calibrated to the stronger wheel session (104638); weaker
intermittent sequences of the kind the paper notes in mouse 102335 could be
missed. The macaque has few units per session (25–64) and a trial-structured task.

Depositing the EBRAINS wheel data on DANDI would be valuable — it is the only
public recording in the condition where the effect is known to exist.

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
  09_validate_windowed.py  # windowed detector validated on both original wheel mice
  10_macaque_control.py    # macaque EC: is the 'oscillation' task-engagement structure?
  11_ca1_rest_000552.py    # CA1 during long TASK-FREE rest (the condition, minus the region)
  12_ca1_rest_figure.py    # the rhythmicity-vs-sequences dissociation + rasters
  13_power_analysis.py     # degrade the positive control to CA1's n and rate: does it still fire?
  scan_001701_regions.py   # metadata scan: which sessions actually have MEC/V1 units
figures/  results/  cache/
replication_notebook.py / .ipynb   # narrative walkthrough
```
