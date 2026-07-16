"""Diagnose the harder positive control (mouse 102335).

The paper computes the neuron sorting on a SUBSET of the session for this mouse
([1,100, 1,400] s) because the sequences are more salient there. We reproduce
that: sort on the subset, then display the whole session with that sorting, for
good-only and good+mua unit sets.
"""
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

fig, axes = plt.subplots(2, 2, figsize=(17, 10))

for row, qual in enumerate([("good",), ("good", "mua")]):
    sess = loaders.load_ebrains("102335", "Wheel-HeadFixed", quality=qual)
    t, z, b = mo.build_activity(sess.spike_times, sess.t0, sess.t1)
    label = "good only" if qual == ("good",) else "good + mua"

    # sort on the paper's subset [1100,1400] s (absolute wheel time; sess.t0 is
    # wheel_start+300). Map to indices in t.
    sub_lo, sub_hi = 1100, 1400
    # times in t are absolute seconds; wheel started at sess.t0-300
    wheel_start = sess.t0 - 300.0
    lo = np.searchsorted(t, wheel_start + sub_lo)
    hi = np.searchsorted(t, wheel_start + sub_hi)
    if hi - lo < 50:               # fallback: middle third
        lo, hi = len(t) // 3, 2 * len(t) // 3

    _, _, theta_sub, *_ = mo._pop_metrics(z[lo:hi])
    order_sub = np.argsort(-theta_sub)
    _, _, theta_full, *_ = mo._pop_metrics(z)
    order_full = np.argsort(-theta_full)

    for col, (order, ttl) in enumerate([(order_full, "sorted on FULL session"),
                                        (order_sub, "sorted on [1100,1400]s subset")]):
        ax = axes[row, col]
        show = min(len(t), int(1500 / mo.BIN_SIZE))
        im = ax.imshow(z[:show, order].T, aspect="auto", cmap="viridis",
                       vmin=-1, vmax=3, origin="lower",
                       extent=[t[0]/60, t[show-1]/60, 0, sess.n_units])
        ax.axvspan((wheel_start+sub_lo)/60, (wheel_start+sub_hi)/60,
                   color="red", alpha=0.12)
        ax.set_title(f"102335 wheel — {label}, {sess.n_units} units — {ttl}", fontsize=10)
        ax.set_xlabel("time (min)"); ax.set_ylabel("cell (sorted)")

fig.suptitle("Mouse 102335 wheel: are periodic sequences visible? (red = paper's sort subset)",
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = f"{ROOT}/figures/diag_102335.png"
fig.savefig(out, dpi=110)
print("->", out)
