"""Cross-dataset summary figure."""

import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ORDER = ["ctrl_wheel", "ctrl_openfield", "d000053_mec", "d001701_mec", "d001701_v1", "d000897_ec"]
LABEL = {
    "ctrl_wheel": "ORIGINAL\nwheel / darkness\n(positive control)",
    "ctrl_openfield": "ORIGINAL\nopen field\n(same units)",
    "d000053_mec": "DANDI 000053\nmouse MEC\nVR track",
    "d001701_mec": "DANDI 001701\nmouse MEC\nX-maze",
    "d001701_v1": "DANDI 001701\nvisual cortex\n(region control)",
    "d000897_ec": "DANDI 000897\nMACAQUE EC\nmental nav.",
}
COL = {"ctrl_wheel": "#2166ac", "ctrl_openfield": "#92c5de", "d000053_mec": "#d6604d",
       "d001701_mec": "#d6604d", "d001701_v1": "#b2b2b2", "d000897_ec": "#d6604d"}

res = {}
for fn in glob.glob(f"{ROOT}/results/*.json"):
    d = json.load(open(fn))
    if d.get("tag") in ORDER:
        res[d["tag"]] = d

tags = [t for t in ORDER if t in res]
x = np.arange(len(tags))
cols = [COL[t] for t in tags]

fig, axes = plt.subplots(1, 3, figsize=(19, 7.2))

# 1. single-cell rhythmicity -- non-specific
ax = axes[0]
vals = [100 * res[t]["frac_rhythmic"] for t in tags]
ax.bar(x, vals, color=cols, edgecolor="k", lw=0.6)
ax.set_xticks(x)
ax.set_xticklabels([LABEL[t] for t in tags], fontsize=7.4, rotation=30, ha="right")
ax.set_ylabel("% of cells with ultraslow rhythmicity")
ax.set_ylim(0, 100)
ax.set_title("(a) Single-cell rhythmicity is NOT specific\n"
             "high everywhere — incl. open field & visual cortex",
             fontsize=10.5, pad=10)
for xi, v in zip(x, vals):
    ax.text(xi, v + 2, f"{v:.0f}%", ha="center", fontsize=8.5)

# 2. population sequence: rotation z  -- the discriminating test
ax = axes[1]
vals = [res[t]["rot_z"] for t in tags]
hatches = ["//" if t == "d001701_v1" else "" for t in tags]
bars = ax.bar(x, vals, color=cols, edgecolor="k", lw=0.6)
for bar, h in zip(bars, hatches):
    bar.set_hatch(h)
ax.axhline(0, color="k", lw=0.8)
ax.axhline(3, color="crimson", ls="--", lw=1.2, label="significance threshold (z=3)")
ax.set_xticks(x)
ax.set_xticklabels([LABEL[t] for t in tags], fontsize=7.4, rotation=30, ha="right")
ax.set_ylabel("population rotation index (z vs circular shift)")
ax.legend(fontsize=8.5, loc="lower left")
ax.set_title("(b) Coherent population sequences\n"
             "in MEC: only in wheel/darkness. Hatched bar = behavioural confound",
             fontsize=10.5, pad=10)
ax.set_ylim(-2.2, 8.2)
ax.annotate("V1 'positive' is the lap cycle,\nnot an intrinsic oscillation\n(see behaviour control)",
            xy=(3.72, 5.4), xytext=(0.45, 6.3), fontsize=8, color="#444",
            arrowprops=dict(arrowstyle="->", color="#444", lw=1))
for xi, v in zip(x, vals):
    ax.text(xi, v + (0.25 if v >= 0 else -0.55), f"{v:+.1f}", ha="center", fontsize=8.5)

# 3. oscillation period of rhythmic cells
ax = axes[2]
periods = [1 / res[t]["median_peak_hz"] if res[t]["median_peak_hz"] else np.nan for t in tags]
ax.bar(x, periods, color=cols, edgecolor="k", lw=0.6)
ax.set_xticks(x)
ax.set_xticklabels([LABEL[t] for t in tags], fontsize=7.4, rotation=30, ha="right")
ax.set_ylabel("median period of rhythmic cells (s)")
ax.axhspan(10, 200, color="grey", alpha=0.12)
ax.set_title("(c) Periods are in the minute-scale range\n"
             "but that alone does not establish the phenomenon",
             fontsize=10.5, pad=10)
for xi, v in zip(x, periods):
    ax.text(xi, v + 3, f"{v:.0f}s", ha="center", fontsize=8.5)

fig.suptitle("Replicating minute-scale oscillatory sequences (Gonzalo Cogno et al. 2024) "
             "in other public entorhinal datasets",
             fontsize=13.5, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = f"{ROOT}/figures/SUMMARY.png"
fig.savefig(out, dpi=125)
print("->", out)

print(f"\n{'session':<16}{'units':>6}{'rhythmic':>10}{'PC12 z':>9}{'rot z':>8}  sequences")
print("-" * 62)
for t in tags:
    r = res[t]
    print(f"{t:<16}{r['n_units']:>6}{100*r['frac_rhythmic']:>9.0f}%{r['pc12_z']:>9.1f}"
          f"{r['rot_z']:>8.1f}  {'YES' if r['sequences'] else 'no'}")
