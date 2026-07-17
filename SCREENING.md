# Screening DANDI for minute-scale oscillatory sequences

Reference appendix for the [README](README.md). What was looked at, what was
rejected, and why. Recorded so nobody repeats the search.

Scope: metadata for **all 875 dandisets**, then streaming
`intervals`/`epochs`/`units`/`electrodes` for every plausible candidate.

## Two findings explain most rejections

**1. The 17.6-minute floor decides almost everything.** The paper estimates PSDs
with Welch, Hamming windows of **17.6 min** (8,192 bins × 129 ms) at 50% overlap,
noting that "a large window was needed to identify oscillation frequencies ≪0.1 Hz".
A single analysis window therefore needs ≥17.6 min of *continuous* recording in one
condition. Public data is almost always chopped finer than that by task structure,
and the datasets that do offer long blocks usually lack sorted units covering them.

**2. The condition is never in the metadata.** Of all 875 dandisets, exactly **3**
mention darkness and none are relevant; **zero** mention "ultraslow", "infraslow" or
"slow oscillation". Where a rest or dark block exists it is discoverable only by
streaming the file, and often is not labelled even there.

## Rejected

| dandiset | what it has | why rejected |
|---|---|---|
| **000582** (Sargolini) | rat MEC, open field | 3–5 units/session |
| **000638** (Shuman) | mouse MEC+HPC, head-fixed VR, up to 158 MEC units, 90–187 min | no `epochs`/`intervals` table at all, so no rest or dark block is identifiable; screen on throughout; 22/33 sessions pilocarpine-epileptic; control-session MEC median only 55 units |
| **000943** (Clark & Nolan) | mouse MEC tetrodes | raw ecephys only, no spike sorting |
| **000409** (IBL Brain Wide Map) | 459 sessions; 19 with ENTm | 639 good MEC units archive-wide (ENTm = 0.49% of 267,264 electrodes); the `spontaneousActivity` block is **exactly 600 s (10.0 min)** by protocol — median and max both 600.1 s, so **no** session can reach 17.6 min |
| **000022** (Allen Visual Coding, FC) | the only genuine **30.04-min contiguous** spontaneous block, in 26/26 sessions (verified: single row, t=4391–6194) | zero ENT units (SUB/ProS/POST only) |
| **000021 / 000713** (Allen) | large Neuropixels panels | zero ENT; spontaneous fragmented, longest single block 5.0 / 4.8 min |
| **000939** (Duszkiewicz/Peyrache) | postsubiculum; `home_cage` epochs 45–100 min; 42–185 units | postsubiculum abuts parasubiculum, the paper's *own* negative control, so a null is expected and uninformative. Separately, 9/31 sessions have broken sleep scoring: bouts starting after the recording ends, and NREM scored during open-field running (one session at 100%) |
| **000056** (Peyrache/Buzsáki) | AD thalamus + postsubiculum, NREM bouts to 28 min | `electrodes.location` is `unknown` for every electrode in every session — AD and PoSub cannot be separated from the file; 16–292 units, only 1 session ≥100 |
| **000987** (Peyrache) | RSC, 5.3 h sessions, scored SWS, UP/DOWN states | 50–71 units, far short of a population test |
| **001695** (Gonzalez/Vöröslakos) | 110–541 units; homecage rest | no entorhinal coverage; `electrodes.location` is a uniform `CA1` placeholder for all 384 channels (true region only in `units.cell_area`) |
| **001703** (lateral EC) | the only LEC dataset on DANDI | **empty — 0 files** |
| **000574 / 000575** (human MTL) | human single units incl. MTL | 28 units/session over 21.5 min; microwire yields cannot support a population test |
| **000230** (Jacobsen et al.) | the only Moser-lab deposit on DANDI | rabies tracing + photostimulation, 9 files, 2 subjects |
| **000059** (Petersen & Buzsáki) | medial septum **cooling**, graded 37.8→19.9 °C — a real causal handle | longest continuous cooled block is **735 s (12.25 min)**; the whole temperature-instrumented recording is 26.1 min. Also hippocampus not MEC, `electrodes.location` `unknown`, trials tile the session, and units tables include noise clusters (one subject lists 1,114 rows of which 40 are `good`) |
| **001375** (Septum GABA DREADDs) | septal chemogenetics | 3-file pilot, 2 mice; VR laps tile the session; raw only |
| **001634** (Bray/Frank, MS optogenetics) | **ten task-free rest epochs of 17.2–33.4 min** ("animal rests in a small empty box", stimulation off) — the best-structured rest data found anywhere; ~1,110 units/session | **units cover only the `lineartrack` epochs; every one of the ten rest epochs has 0 units.** Units are Spyglass-style, one per file across 434 files, so the main session file reports `units: False`. Using the rest blocks means sorting 242 GB/session |
| **Frank lab** (000115, 001836, 000629, 001280, 000410, 001566, 000065) | cleanly labelled rest epochs (000115 `02_r1` = 90.3 min; 001836 `01_Seq2Sleep1` = 31.9 min) | **no `units` table in any file** — raw e-series + behavior only |

## Two bugs in this screen

Both would silently mislead anyone repeating it.

1. **Abbreviations defeat keyword search.** 001634 never spells out "medial septum";
   its description says only "PV+ neurons within the **MS**". A regex on the
   spelled-out form misses the best septum dataset on the archive. The inverse also
   bites: matching `mec` as a substring hits "**mec**hanical" and "**mec**hanism".
2. **Dandiset metadata systematically understates content.** 000690's description
   never mentions entorhinal cortex, yet it holds the archive's best MEC coverage and
   was missed by the original triage entirely. 000552's epochs table carries only
   start/stop times, so its 160-min task-free blocks are invisible unless trial
   coverage is checked empirically. Any screen that reads descriptions instead of
   streaming files is wrong in both directions.

## The pattern

Derived data covers what the original paper analyzed, not what a reuser needs. 001634
sorted the running epochs and left ten ideal rest epochs unsorted; the Frank lab sets
deposit no units at all; 000053's lightweight companion files exist only for the
non-darkness sessions. The raw data is there in every case.

For an archive this argues for sorting whole sessions rather than analyzed intervals,
and for making the epoch coverage of derived data discoverable without streaming 434
files to find out.

## Connected regions: candidates and their data

The paper's model splits the phenomenon: the **oscillation** is proposed to be global
and ascending, the **sequence** MEC-local. That predicts a dissociation — other
targets of the same ascending drive should show the oscillation without the sequence.
Ultraslow oscillations are, by the paper's own citations, already reported in
hippocampus (Penttonen 1999), basal ganglia (Ruskin 1999), human entorhinal cortex
(Aghajan, Kreiman & Fried 2023), cerebral cortex (Aladjalova 1957), monkey visual
cortex (Leopold 2003), thalamocortical sleep (Lecci 2017), pupil (Blasiak 2013) and
vasculature (Drew 2020).

| region | why a candidate | public data |
|---|---|---|
| **Lateral entorhinal cortex** | the strongest untested candidate: adjacent, reciprocally connected, and the paper itself invokes "slowly drifting neural population activity in lateral entorhinal cortex". Never recorded in the paper | **none** — 001703 is empty |
| **Presubiculum** | HD input to MEC layer III; PaS was the paper's control, PreS untested | 000939 is postsubiculum, adjacent to the paper's own PaS control |
| **Medial septum / diagonal band** | theta pacemaker and major MEC input — the way to ask whether the rhythm is *imposed* | 000059, 001375, 001634 — none usable as deposited |
| **Brainstem neuromodulatory** (LC, DR, VTA, PPT/LDT) | the paper's own proposed source | no in vivo population recordings on DANDI |
| **AD thalamus / HD ring** | a ring attractor that sweeps during sleep; HD drift is a random walk, *not* periodic, so it would dissociate sequence from oscillation | 000056 (regions inseparable, units too few) |
| **Human EC** | Aghajan et al. report the same phenomenon in human EC — though their claim is single-cell periodicity, which this repo shows is non-specific | 000574/000575: 28 units/session |
