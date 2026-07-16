"""Choose the detector using ground truth, not intuition.

The same 400+ MEC units from the original paper were recorded in two epochs:
    Wheel-HeadFixed : rotating wheel in darkness -> effect KNOWN PRESENT
    OpenField       : free foraging              -> effect UNKNOWN

CORRECTION: the OpenField epoch is NOT a validated negative control. The paper's
Methods state these trials "were not used in the present study", so nothing is
known about the effect there. The detector is therefore selected on the WHEEL
column alone -- the only condition with a known-present effect. The open-field
column is reported as an observation, not as a specificity criterion.

(Specificity evidence comes from the paper itself, Fig. 5: a fraction of cells in
parasubiculum and visual cortex were also "ultraslow and periodic", yet in neither
region were those oscillations organized into sequences. Single-cell rhythmicity
is non-specific in the paper's own data; only the population test discriminates.)

We compare three candidate signals to run the spectrum on:

    A. binarized activity        (paper's literal activity matrix)
    B. smoothed z-scored rate    (no binarization)
    C. z-scored spike-train autocorrelogram  (reference implementation's choice)

and fix the permutation resolution floor by pooling the null across cells:
under H0 the normalized statistic is exchangeable across cells, so 200
surrogates x N cells gives ~1e-5 p-value resolution instead of 1/201.
"""

import os
import sys

import numpy as np
from scipy import signal, stats
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(__file__))
import loaders
import mec_oscillation as mo

N_SURR = 100
BAND = mo.BAND
FS = 1.0 / mo.BIN_SIZE


def acg_zscored(counts, window=560.0, bin_size=mo.BIN_SIZE):
    """Z-scored spike-train autocorrelogram, as in the reference notebook.

    The pairwise-lag correlogram at 120 ms bins is exactly the autocorrelation
    of the 120 ms binned count vector, so we compute it by FFT instead of
    enumerating ~17M spike pairs per cell.
    """
    nlag = int(window / bin_size)
    x = counts - counts.mean()
    full = signal.correlate(x, x, mode="full", method="fft")
    mid = len(x) - 1
    a = full[max(0, mid - nlag): mid + nlag + 1]
    sd = a.std()
    return (a - a.mean()) / (sd if sd > 0 else 1.0)


def binned_counts(spike_times, t0, t1, bin_size=mo.BIN_SIZE):
    edges = np.arange(t0, t1 + bin_size, bin_size)
    return np.array([np.histogram(st, bins=edges)[0].astype(float)
                     for st in spike_times])


def spectra(spike_times, t0, t1, variant):
    """Return (freqs_in_band, psd_in_band[n_cells, n_f]) for a detector variant."""
    if variant == "C":
        cnt = binned_counts(spike_times, t0, t1)
        psds = []
        for c in cnt:
            a = acg_zscored(c)
            f, p = signal.welch(a, FS, window="hamming",
                                nperseg=min(len(a), 8196), detrend="constant")
            psds.append(p)
        f = np.asarray(f)
        psds = np.array(psds)
    else:
        _, z, b = mo.build_activity(spike_times, t0, t1)
        x = b if variant == "A" else z
        psds = []
        for i in range(x.shape[1]):
            f, p = mo.compute_psd(x[:, i], FS)
            psds.append(p)
        f = np.asarray(f)
        psds = np.array(psds)

    fb = (f >= BAND[0]) & (f <= BAND[1])
    return f[fb], psds[:, fb]


def run_variant(sess, variant, seed=0):
    rng = np.random.default_rng(seed)
    f_in, psd_obs = spectra(sess.spike_times, sess.t0, sess.t1, variant)
    n_cells = len(sess.spike_times)

    psd_sur = np.zeros((N_SURR, n_cells, len(f_in)))
    for j in tqdm(range(N_SURR), desc=f"  {variant} surrogates", leave=False):
        sur = [mo.isi_shuffle(st, sess.t0, sess.t1, rng) for st in sess.spike_times]
        _, psd_sur[j] = spectra(sur, sess.t0, sess.t1, variant)

    eps = 1e-30
    tot = psd_sur.sum(axis=0)
    mean_sur = tot / N_SURR
    stat_obs = (psd_obs / np.maximum(mean_sur, eps)).max(axis=1)

    # pooled null across cells (exchangeable under H0) -> fine p resolution
    null = []
    for j in range(N_SURR):
        loo = (tot - psd_sur[j]) / (N_SURR - 1)
        null.append((psd_sur[j] / np.maximum(loo, eps)).max(axis=1))
    null = np.concatenate(null)

    p = np.array([(1 + (null >= s).sum()) / (1 + len(null)) for s in stat_obs])
    sig = mo.benjamini_hochberg(p, 0.05)
    return stat_obs, p, sig, null


if __name__ == "__main__":
    print("Loading the two ground-truth epochs (same units, same probe)...")
    wheel = loaders.load_ebrains("104638", "Wheel-HeadFixed")
    field = loaders.load_ebrains("104638", "OpenField")
    print(f"  wheel     : {wheel.n_units} units, {wheel.duration/60:.1f} min")
    print(f"  open field: {field.n_units} units, {field.duration/60:.1f} min")

    print(f"\n{'variant':<48} {'wheel':>16} {'open field':>16}")
    print("-" * 84)
    names = {
        "A": "A. binarized activity (paper literal)",
        "B": "B. smoothed z-scored rate",
        "C": "C. z-scored spike-train autocorrelogram",
    }
    for v in ("A", "B", "C"):
        sw, pw, gw, nullw = run_variant(wheel, v)
        sf, pf, gf, nullf = run_variant(field, v)
        print(f"{names[v]:<48} {100*gw.mean():>14.0f}% {100*gf.mean():>15.0f}%")
        print(f"{'   median excess ratio':<48} {np.median(sw):>14.2f}x {np.median(sf):>15.2f}x")
        print(f"{'   null 95th pct':<48} {np.percentile(nullw,95):>14.2f}x "
              f"{np.percentile(nullf,95):>15.2f}x")
    print("\nSelect on the WHEEL column: the detector must recover the effect where it is")
    print("known to exist (paper: 78% of Neuropixels units, 683/879). The open-field column")
    print("is an untested condition, NOT a specificity criterion -- see module docstring.")
