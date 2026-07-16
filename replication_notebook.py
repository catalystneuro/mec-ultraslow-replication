# %% [markdown]
# # Minute-scale oscillatory sequences in entorhinal cortex — a replication
#
# Testing **Gonzalo Cogno et al. (2024)**, *Minute-scale oscillatory sequences in
# medial entorhinal cortex*, Nature 625:338–344, on public datasets other than
# the original. Prompted by [dandi/helpdesk#156](https://github.com/dandi/helpdesk/discussions/156).
#
# **The claim.** MEC activity organizes into ultraslow oscillations (periods of tens
# of seconds to minutes), during which cells fire in periodic sequences.
#
# **The condition.** The effect was observed in a *sensory-minimized* setting: mice
# head-fixed on a rotating wheel **in darkness**, no rewards, no change in location or
# heading.
#
# **Darkness is where the effect was observed — NOT a condition shown to be required.**
# That was the only condition the paper ever ran; it contains no light or navigation
# comparison. The authors say so directly: *"It remains an open question whether the
# ultraslow oscillatory sequences are present across a broader spectrum of behaviours,
# including sleep and free exploration, and in the presence of salient visual feedback.
# If so, it is possible that the sequences reset in the presence of strong landmarks or
# sensory stimulation."*
#
# So testing DANDI's navigation datasets is a partial **answer to the paper's own open
# question**, not a defective replication.

# %%
import sys
sys.path.insert(0, "scripts")

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

import loaders
import mec_oscillation as mo

# %% [markdown]
# ## 1. A positive control is mandatory
#
# A null result is worthless unless the detector demonstrably fires when the effect
# is present. The original Neuropixels data is public on EBRAINS
# ([10.25493/SKKX-4W3](https://doi.org/10.25493/SKKX-4W3)) and — conveniently — the
# NWB file holds **two epochs from the same 400+ MEC units**:
#
# * `Wheel-HeadFixed` — rotating wheel in darkness → effect **known present**
# * `OpenField` — free foraging → the Methods state these trials "were not used in the
#   present study", so the effect there is **unknown**. This is an untested condition,
#   *not* a validated negative control — our result on it is a new observation, not a
#   check on the pipeline.
#
# The wheel epoch is the real control. We stream it (no bulk download).

# %%
wheel = loaders.load_ebrains("104638", "Wheel-HeadFixed")
print(f"{wheel.n_units} MEC units, {wheel.duration/60:.1f} min, {wheel.region}")

t, z, b = mo.build_activity(wheel.spike_times, wheel.t0, wheel.t1)
print(f"activity matrix: {z.shape} (time x cells) at {mo.BIN_SIZE*1000:.0f} ms bins, "
      f"Gaussian sigma = {mo.SIGMA} s")

# %% [markdown]
# ### Two statistical traps (both caught by the positive control)
#
# Each of these produced a confident, clean, and completely **wrong null** on data
# where the effect is known to be present.
#
# **Trap 1 — band power cannot detect the effect.** The σ = 5 s Gaussian kernel *is*
# a low-pass filter, so ~75% of the power lands in the ultraslow band for the real
# data *and* for the shuffled surrogates alike (0.74 vs 0.77). A power-in-band
# statistic is blind by construction. The signature is a **narrowband peak above the
# broadband 1/f background**, so the statistic must whiten that background:
# `max over band of PSD_data(f) / mean_surrogate_PSD(f)`.
#
# **Trap 2 — the permutation resolution floor.** With 200 surrogates the smallest
# attainable p-value is 1/201 ≈ 0.005, but Benjamini–Hochberg across 404 cells needs
# p < 0.05/404 ≈ 1.2e-4. *No cell can ever reach significance* regardless of effect
# size — the original data returned **0/404** with a max excess-power ratio of 47×.
# Fixed by pooling the null across cells (the normalized statistic is exchangeable
# under H0), which lifts resolution to ~1e-5.

# %%
sc = mo.single_cell_test(wheel.spike_times, wheel.t0, wheel.t1,
                         n_surrogates=100, progress=tqdm)
print(f"ultraslow-rhythmic cells: {sc['sig'].sum()}/{wheel.n_units} "
      f"= {100*sc['sig'].mean():.0f}%   (paper: ~91% by a different criterion)")
fp = sc["f_peak"][sc["sig"]]
print(f"median period: {1/np.median(fp):.0f} s  -> 'tens of seconds to minutes' ✓")

# %% [markdown]
# We also had to pick *which signal* to run the spectrum on empirically, using the
# wheel-vs-open-field ground truth (`03_detector_bakeoff.py`):
#
# | detector | wheel (effect present) | open field (same cells) |
# |---|---|---|
# | binarized activity (paper's literal matrix) | 2% | 0% |
# | smoothed z-scored rate | 0% | 0% |
# | **z-scored autocorrelogram spectrum** | **64%** | **73%** |
#
# The detector is selected on the **wheel** column (the known-present condition):
# only the autocorrelogram spectrum recovers the effect. The oscillation's phase is
# non-stationary, so Welch on the raw trace averages segments and smears the peak
# away, while the autocorrelogram averages over the whole session and survives it.
# (Compare the paper's own Neuropixels figure: 78% of units, 683/879.)
#
# The second column is informative but is *not* a specificity check (open field is an
# untested condition). The real specificity evidence comes from the paper itself,
# Fig. 5: *"while the calcium activity of a fraction of cells in both brain areas
# [parasubiculum and visual cortex] was ultraslow and periodic, in neither brain region
# were these oscillations organized into oscillatory sequences."*
#
# **A single-cell ultraslow spectral peak is therefore not specific evidence for this
# phenomenon — in the paper's data either.** Any replication that stops there reports a
# false positive. Only the *population sequence* test discriminates.

# %% [markdown]
# ## 2. The real test: coherent population sequences
#
# The paper's actual claim is that cells ride a **common** oscillation at **staggered
# phases**, forming a repeating sequence. The right null therefore preserves each
# cell's own slow modulation and destroys only the coordination *between* cells: an
# independent **circular shift** of each cell's trace. A circular shift leaves every
# cell's power spectrum and autocorrelation exactly intact.

# %%
pop = mo.population_sequence_test(z[::8], n_surrogates=200, progress=tqdm)
print(f"PC1+PC2 variance : {100*pop['pc12_var']:.1f}%  vs {100*pop['null_pc12_mean']:.1f}% shuffled")
print(f"rotation index   : {pop['rotation']:.3f} vs {pop['null_rot_mean']:.3f} shuffled "
      f"(z = {pop['rot_z']:+.1f})")

# %% [markdown]
# Sorting cells by the angle of their (PC1, PC2) loading vector — exactly as the
# paper does — exposes the periodic sequences as repeating diagonal bands.

# %%
order = np.argsort(-pop["theta"])
fig, ax = plt.subplots(figsize=(15, 5))
show = int(1200 / mo.BIN_SIZE)
im = ax.imshow(z[:show, order].T, aspect="auto", cmap="viridis", vmin=-1, vmax=3,
               origin="lower", extent=[t[0]/60, t[show-1]/60, 0, wheel.n_units])
ax.set_xlabel("time (min)"); ax.set_ylabel("cell (sorted by PCA angle θ)")
ax.set_title("Positive control: periodic sequences in MEC on the wheel, in darkness")
plt.colorbar(im, ax=ax, label="z-scored rate")
plt.show()

# %% [markdown]
# ## 3. Apply the validated pipeline to the other datasets
#
# Screened from DANDI: **000053** (mouse MEC, Neuropixels, VR track), **001701**
# (mouse MEC, X-maze — and *visual cortex on the same probe*, reproducing the
# paper's own regional control), **000897** (**macaque** EC, mental navigation).
# Rejected: 000582 (3–5 units/session), 000638 (too few MEC units), 000943 (no
# spike sorting).
#
# For robustness we run EVERY session (`07_run_all_sessions.py`, resumable) and
# report the paper's own quantity — the FRACTION OF SESSIONS with sequences
# (`08_session_summary.py`, `figures/SUMMARY_all_sessions.png`):
#
# | group | condition | sessions with sequences |
# |---|---|---|
# | **EBRAINS wheel / darkness** | positive control | **2/2 = 100%** |
# | EBRAINS open field (same units) | untested in paper | 0/2 |
# | **DANDI 000053** mouse MEC | VR track | **0/20 = 0%** |
# | **DANDI 001701** mouse MEC | X-maze | **3/114 = 3%** |
# | DANDI 001701 visual cortex | region control | 9/110 = 8% *(confound)* |
# | DANDI 000897 macaque EC | mental navigation | 8/15 = 53% *(confound)* |
#
# The detector fires on BOTH original wheel sessions and is deliberately not tuned
# to force borderline cases. Single-session `04_analyze.py` still exists for the
# detailed per-session figures.

# %% [markdown]
# ## 4. Behaviour/task structure manufactures the signature — in two species
#
# The two apparent positives outside the wheel are both behavioural artifacts:
#
# - **Mouse visual cortex** (9/110, higher than MEC's 3/114) — the region the paper
#   found negative. The animal runs laps, the scene repeats, V1 tiles the lap cycle.
#   `05_behavior_control.py`: position, MEC-PC1 and V1-PC1 all peak at 0.007–0.010 Hz.
# - **Macaque EC** (8/15 — at first glance a striking cross-species replication).
#   `10_macaque_control.py`: population PC1 tracks trial-onset density (drops to zero
#   during task pauses, resumes with the trials); its spectrum overlaps the
#   trial-onset spectrum across the band. It is task-engagement structure over a
#   multi-hour session, not an intrinsic rhythm.

# %% [markdown]
# ## Conclusion
#
# The pipeline **reproduces the original finding** on the wheel data (both Neuropixels
# mice), and finds **no minute-scale oscillatory sequences in mouse MEC during
# navigation** across 134 sessions (0/20 VR + 3/114 X-maze). The apparent positives in
# mouse V1 (9/110) and macaque EC (8/15) are behavioural/task confounds, each nailed by
# a dedicated control.
#
# This is consistent with the authors' own speculation that "sequences reset in the
# presence of strong landmarks or sensory stimulation", and is a partial answer to the
# open question they pose. It does **not** establish that darkness is *necessary* — the
# paper never tested another condition, and our nulls cannot separate "sensory drive
# abolishes the sequences" from "insufficient power to detect them".
#
# **The biggest caveat:** the paper found oscillatory sequences in only **15 of 27** MEC
# wheel sessions (Fig. 5g). The effect is highly session-variable even under ideal
# conditions, and we ran **one session per dataset**. The clear next step is to run all
# sessions per dandiset (001701 has 218) and report the *fraction of sessions* with
# sequences — the quantity the paper actually reports.
#
# Depositing the EBRAINS wheel data on DANDI would still be valuable: it is the only
# public recording in the condition where the effect is known to exist.
