"""CA1 task-free rest (000552) vs the paper's own wheel data.

Makes the dissociation explicit: CA1 at rest is the MOST single-cell rhythmic
group in the whole panel -- more so than the positive control where the effect is
known present -- and yet has no population sequences at all. That is the cleanest
demonstration in this repo that single-cell ultraslow rhythmicity is not evidence
for the paper's claim, and that only the population test discriminates.
"""

import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import loaders
import mec_oscillation as mo

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DECIM = 8

# dy = label offset (points), to keep crowded labels legible
GROUPS = [
    ("ebrains", "Wheel-HeadFixed", "ORIGINAL wheel\n(darkness)", "#2166ac", 15),
    ("ebrains", "OpenField", "ORIGINAL\nopen field", "#92c5de", 15),
    ("000053", "MEC", "000053 MEC\nVR track", "#d6604d", 15),
    ("001701", "MEC", "001701 MEC\nX-maze", "#d6604d", -38),
    ("001701", "V1", "001701 V1", "#b2b2b2", 15),
    ("000897", "EC", "000897 macaque EC", "#b2b2b2", 15),
    ("000552", "CA1-Awake", "000552 CA1\ntask-free rest", "#5aae61", 15),
]


def collect():
    recs = {}
    for fn in glob.glob(f"{ROOT}/results/sessions/*.json"):
        d = json.load(open(fn))
        if d.get("status") != "ok":
            continue
        recs.setdefault((d["dandiset"], d["region_select"]), []).append(d)
    return recs


def best_window(z, win_s=500.0, step_s=100.0):
    """Start time (s) of the window with the highest rotation index."""
    bin_s = mo.BIN_SIZE * DECIM
    w, step = int(win_s / bin_s), int(step_s / bin_s)
    best = (0.0, 0.0)
    for i0 in range(0, max(1, z.shape[0] - w), step):
        _, rot, *_ = mo._pop_metrics(z[i0:i0 + w])
        if rot > best[1]:
            best = (i0 * bin_s, rot)
    return best


def raster_panel(ax, z, title, win_s=500.0, start_s=0.0):
    """Sort cells on a WINDOW and display that same window.

    The paper sorts its own Neuropixels rasters on a subset ([1200,1700] s for
    mouse 104638) rather than the whole session, because the sequences are
    intermittent and a whole-session sort smears them. We do the same, and apply
    the identical procedure to both panels so the comparison is fair.
    """
    bin_s = mo.BIN_SIZE * DECIM
    i0 = int(start_s / bin_s)
    i1 = min(i0 + int(win_s / bin_s), z.shape[0])
    seg = z[i0:i1]
    x = seg - seg.mean(0, keepdims=True)
    _, _, theta, _, _, _ = mo._pop_metrics(x)
    order = np.argsort(-theta)
    b = (x > 1.0).astype(float)[:, order].T
    ax.imshow(b, aspect="auto", cmap="Greys", interpolation="nearest",
              extent=[start_s, start_s + b.shape[1] * bin_s, b.shape[0], 0])
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Neuron (sorted on this window)")


def main():
    recs = collect()

    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1])

    # ---- (a) the dissociation
    ax = fig.add_subplot(gs[0, :])
    for dsid, reg, label, col, dy in GROUPS:
        rs = recs.get((dsid, reg), [])
        if not rs:
            continue
        rhy = 100 * np.mean([r["frac_rhythmic"] for r in rs])
        seq = 100 * np.mean([r["sequences"] for r in rs])
        ax.scatter(rhy, seq, s=190, color=col, edgecolor="k", zorder=3)
        ax.annotate(f"{label}\n({sum(r['sequences'] for r in rs)}/{len(rs)})",
                    (rhy, seq), textcoords="offset points", xytext=(0, dy),
                    ha="center", fontsize=8.5)
    ax.set_xlabel("single-cell ultraslow rhythmicity (% of cells)  —  NON-SPECIFIC")
    ax.set_ylabel("% sessions with\nperiodic sequences")
    ax.set_xlim(20, 100)
    ax.set_ylim(-12, 118)
    ax.axhline(0, color="grey", lw=0.6, ls=":")
    ax.set_title("(a) Single-cell rhythmicity does not predict population sequences.\n"
                 "CA1 at rest is the MOST rhythmic group in the panel — more than the "
                 "positive control — yet has zero sequences.", fontsize=11, pad=10)

    # ---- (b, c) rasters: positive control vs CA1 rest, each at ITS OWN best
    # window (highest rotation index). Showing each dataset at its best case is
    # the fair comparison: it gives CA1 every chance to reveal a sequence.
    try:
        s = loaders.load_ebrains("104638", "Wheel-HeadFixed")
        z = mo.build_activity(s.spike_times, s.t0, s.t1)[1][::DECIM]
        st, rot = best_window(z)
        raster_panel(fig.add_subplot(gs[1, 0]), z,
                     f"(b) POSITIVE CONTROL — original MEC, wheel in darkness\n"
                     f"{s.n_units} units — best window (rotation={rot:.2f}); "
                     f"sequences present, p=1e-6", start_s=st)
    except Exception as e:
        print("ebrains raster unavailable:", e)

    s2 = loaders.load_000552(
        "sub-e15-13f1/sub-e15-13f1_ses-e15-13f1-220119_behavior+ecephys.nwb")
    z2 = mo.build_activity(s2.spike_times, s2.t0, s2.t1)[1][::DECIM]
    st2, rot2 = best_window(z2)
    raster_panel(fig.add_subplot(gs[1, 1]), z2,
                 f"(c) CA1, task-free awake rest (000552)\n"
                 f"{s2.n_units} units, 86 min, 90% rhythmic — best window "
                 f"(rotation={rot2:.2f}); no sequences", start_s=st2)

    fig.suptitle("Hippocampal CA1 during long task-free rest: rhythmic but not sequential",
                 fontsize=13.5, y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"{ROOT}/figures/ca1_rest_000552.png"
    fig.savefig(out, dpi=125)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
