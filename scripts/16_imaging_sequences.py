"""The population sequence test on the paper's PRIMARY data: two-photon imaging.

This closes the biggest remaining gap in the panel. Every sequence result here and in
the paper's own validation rests on **two Neuropixels mice**. But the paper's headline
numbers -- 91% of cells oscillatory, ~94% phase-locked, 44% at 0.006-0.008 Hz, the
ring manifold, 15/27 oscillatory sessions -- all come from **two-photon imaging of
6,231 cells across 15 sessions**. That is the actual claim, and no population sequence
test had touched it.

The single-cell half of the imaging data already replicates independently: see
https://github.com/rly/replicate-gonzalo-cogno-2023, which reproduces Figure 1 on
this exact session (mouse 60584, session 7 -- the paper's own Fig. 1b/2a-e example),
recovering the reported ~0.0066 Hz peak and its harmonics. But single-cell rhythmicity
is precisely what this repo shows is NON-SPECIFIC (83% of CA1 cells at rest are
rhythmic with zero sequences; in MEC it is pupil-coherent in 8/11 sessions). So
reproducing Fig. 1 does not test the distinctive claim. The population test does.

Data: EBRAINS doi:10.25493/SKKX-4W3, data/calcium/60584/2019-01-29/MUnit_0/suite2p.

Preprocessing follows the paper's "Binary deconvolved calcium activity" methods:
Suite2p deconvolved (spks), manual curation (iscell), downsample x4 by mean to
129 ms / 7.73 Hz, binarize at mean + 1.5 SD.

  Cell count note: the paper reports 484 cells for this session; iscell gives 522.
  The SNR>4 criterion is non-binding on this session -- reimplementing it on Fcorr
  (F - 0.7*Fneu) gives a median SNR of 21.5 and NO cell in 0 < SNR <= 4, so it
  excludes nobody with a valid noise estimate. The 7% gap cannot matter for a
  population test in which the paper says ~94% of cells participate.

DETECTOR INPUT -- and a trap. The primary variant is the paper's own activity matrix
(binarize the deconvolved trace at mean + 1.5 SD, no extra smoothing). The ephys
pipeline's convention (Gaussian sigma = 5 s, then z-score) is ALSO run, and it must
not be used here: deconvolved calcium is already temporally smooth, so smoothing again
oversmooths, and the PC1/PC2 trajectory then advances consistently for the circular-
shift SURROGATES too (observed rotation 0.853 vs a null of ~0.83, z=+1.2). The
statistic saturates and loses all discriminative power, returning a confident false
negative. The sigma = 5 s kernel exists in the ephys path to mimic calcium dynamics
that imaging data already has -- applying it to calcium double-counts it. This is why
the paper binarizes imaging data instead.

usage: python3 scripts/16_imaging_sequences.py
"""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import requests
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(__file__))
import mec_oscillation as mo

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = f"{ROOT}/cache/calcium_60584_2019-01-29_MUnit_0"
EB = ("https://rgw.cscs.ch/ebrains:d-a25f722e-f9b1-4f32-80bc-37206c2d1536"
      "/data/calcium/60584/2019-01-29/MUnit_0/suite2p/plane0")

FS_RAW = 30.95           # Hz, paper
DOWNSAMPLE = 4           # -> 7.73 Hz, 129 ms bins
BIN_S = DOWNSAMPLE / FS_RAW
DECIM = 8                # -> ~1.03 s, matching the ephys pipeline's 0.96 s
N_SHIFT = 200


def fetch(name):
    os.makedirs(CACHE, exist_ok=True)
    p = f"{CACHE}/{name}"
    if not os.path.exists(p):
        print(f"  downloading {name} ...", flush=True)
        r = requests.get(f"{EB}/{name}", timeout=1800)
        r.raise_for_status()
        open(p, "wb").write(r.content)
    return np.load(p)


def preprocess():
    spks = fetch("spks.npy")
    iscell = fetch("iscell.npy")
    keep = iscell[:, 0] == 1
    d = spks[keep].astype(float)
    n_bins = d.shape[1] // DOWNSAMPLE
    d = d[:, : n_bins * DOWNSAMPLE]
    d = d.reshape(d.shape[0], n_bins, DOWNSAMPLE).mean(axis=2)   # paper: mean over ~129 ms
    print(f"  {spks.shape[0]} ROIs -> {keep.sum()} curated cells; "
          f"{n_bins} bins @ {BIN_S*1000:.0f} ms = {n_bins*BIN_S/60:.1f} min")
    return d


def main():
    print("Loading the paper's own Fig. 1b/2a-e session (mouse 60584, session 7)")
    dec = preprocess()

    # (a) paper's activity matrix: binarize the deconvolved trace at mean + 1.5 SD
    thr = dec.mean(1, keepdims=True) + 1.5 * dec.std(1, keepdims=True)
    binary = (dec > thr).astype(float).T                    # (time, cells)

    # (b) detector input, mirroring build_activity for ephys
    sm = gaussian_filter1d(dec, sigma=5.0 / BIN_S, axis=1)
    z = ((sm - sm.mean(1, keepdims=True)) /
         np.where(sm.std(1, keepdims=True) == 0, 1, sm.std(1, keepdims=True))).T

    out = {}
    for label, mat, bs in (("paper's binary matrix", binary[::DECIM], BIN_S * DECIM),
                           ("ephys smoothing convention (INVALID here - oversmooths)",
                            z[::DECIM], BIN_S * DECIM)):
        print(f"\n--- {label}: {mat.shape[0]} bins x {mat.shape[1]} cells")
        pop = mo.population_sequence_test(mat, n_surrogates=N_SHIFT)
        win = mo.windowed_sequence_test(mat, bs, n_surrogates=N_SHIFT)
        print(f"    whole-session: pc12_var={pop['pc12_var']:.3f} (z={pop['pc12_z']:+.1f}) "
              f"rotation={pop['rotation']:.3f} (z={pop['rot_z']:+.1f})")
        print(f"    WINDOWED     : {win['n_sig_windows']}/{win['n_windows']} windows "
              f"significant, oscillation score={win['frac_sig_windows']:.2f}, "
              f"p_session={win['p_session']:.2g} -> sequences={win['sequences']}")
        out[label] = dict(pc12_var=pop["pc12_var"], pc12_z=pop["pc12_z"],
                          rotation=pop["rotation"], rot_z=pop["rot_z"],
                          n_windows=win["n_windows"], n_sig=win["n_sig_windows"],
                          oscillation_score=win["frac_sig_windows"],
                          p_session=win["p_session"], sequences=bool(win["sequences"]))
    out["n_cells"] = int(binary.shape[1])
    out["duration_min"] = float(binary.shape[0] * BIN_S / 60)
    json.dump(out, open(f"{ROOT}/results/imaging_sequences.json", "w"), indent=2)

    # ---------------- figure: the paper's Fig. 2b/2c on its own data ----------------
    # Raster sorting: the paper's binary activity matrix (Fig. 2b).
    # Manifold: the paper convolves each row with a Gaussian before PCA for the
    # manifold specifically ("Manifold visualization for MEC sessions"), so the ring
    # is visualized on the smoothed matrix (Fig. 2c). Smoothing is right for the
    # PICTURE and wrong for the STATISTIC -- see the trap in the module docstring.
    _, _, theta, _, _, _ = mo._pop_metrics(binary[::DECIM])
    order = np.argsort(-theta)
    _, _, _, _, scores, ve = mo._pop_metrics(z[::DECIM])
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))

    show = binary[:, order][: int(3000 / BIN_S)].T
    ax[0].imshow(show, aspect="auto", cmap="Greys", interpolation="nearest",
                 extent=[0, show.shape[1] * BIN_S, show.shape[0], 0])
    ax[0].set(title="(a) PCA-sorted raster — the paper's Fig. 2b\n"
                    f"{binary.shape[1]} cells, mouse 60584 session 7",
              xlabel="Time (s)", ylabel="Neuron (PCA-sorted)")

    ax[1].plot(scores[:, 0], scores[:, 1], lw=0.3, alpha=0.6, color="C0")
    ax[1].set(title=f"(b) PC1–PC2 manifold — the paper's Fig. 2c (smoothed, as they do)\n"
                    f"PC1={ve[0]*100:.1f}% PC2={ve[1]*100:.1f}%  ratio={ve[0]/ve[1]:.2f}"
                    f"  (ring => ~1.0)",
              xlabel="PC1", ylabel="PC2")
    ax[1].set_aspect("equal")

    r = out["paper's binary matrix"]
    ax[2].axis("off")
    ax[2].text(0.02, 0.95, "(c) Population sequence test on the PRIMARY data",
               fontsize=12, weight="bold", va="top")
    bad = out["ephys smoothing convention (INVALID here - oversmooths)"]
    txt = (f"cells: {out['n_cells']}   duration: {out['duration_min']:.0f} min\n"
           f"preprocessing: the paper's own\n(deconvolved, binarized at mean+1.5 SD)\n\n"
           f"windowed sequence test:\n"
           f"   significant windows : {r['n_sig']}/{r['n_windows']}\n"
           f"   oscillation score   : {r['oscillation_score']:.2f}\n"
           f"   p_session           : {r['p_session']:.2g}\n"
           f"   SEQUENCES           : {'YES' if r['sequences'] else 'NO'}\n\n"
           f"the paper's own Neuropixels\nwheel sessions, for reference:\n"
           f"   104638 osc-score 0.54, p=1e-6\n"
           f"   102335 osc-score 0.23, p=0.025\n\n"
           f"TRAP: smoothing is right for the\nmanifold PICTURE (panel b) and\n"
           f"wrong for the STATISTIC: with\nsigma=5 s the null saturates ->\n"
           f"   rotation {bad['rotation']:.2f} but z={bad['rot_z']:+.1f}\n"
           f"   p={bad['p_session']:.2g}: a false negative")
    ax[2].text(0.02, 0.82, txt, fontsize=10.5, va="top", family="monospace")

    fig.suptitle("Minute-scale oscillatory SEQUENCES in the paper's primary two-photon data",
                 fontsize=13.5, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    p = f"{ROOT}/figures/imaging_sequences.png"
    fig.savefig(p, dpi=125)
    print(f"\n-> {p}")


if __name__ == "__main__":
    main()
