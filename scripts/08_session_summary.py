"""Aggregate all per-session results -> fraction of SESSIONS with sequences.

This is the quantity the paper reports (Fig. 5g): MEC 15/27 sessions oscillatory,
PaS 0/25, VIS 0/19. A per-session fraction is robust to the effect's strong
session-to-session variability in a way a single session is not.
"""

import glob
import json
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

GROUPS = [
    ("ebrains", "Wheel-HeadFixed", "ORIGINAL\nwheel / darkness\n(effect known present)", "#2166ac"),
    ("ebrains", "OpenField", "ORIGINAL\nopen field\n(unanalyzed in paper)", "#92c5de"),
    ("000053", "MEC", "DANDI 000053\nmouse MEC\nVR track", "#d6604d"),
    ("001701", "MEC", "DANDI 001701\nmouse MEC\nX-maze", "#d6604d"),
    ("001701", "V1", "DANDI 001701\nvisual cortex\n(region control)", "#b2b2b2"),
    ("000897", "EC", "DANDI 000897\nMACAQUE EC\nmental nav.", "#d6604d"),
    ("000552", "CA1-Awake", "DANDI 000552\nmouse CA1\ntask-free rest\n(region control)", "#b2b2b2"),
]

recs = defaultdict(list)
status = defaultdict(lambda: defaultdict(int))
for fn in glob.glob(f"{ROOT}/results/sessions/*.json"):
    d = json.load(open(fn))
    g = (d.get("dandiset"), d.get("region_select"))
    status[g][d.get("status")] += 1
    if d.get("status") == "ok":
        recs[g].append(d)

print(f"{'group':<34}{'ok':>5}{'skip':>6}{'err':>5}   sessions with sequences")
print("-" * 84)

rows = []
for dsid, reg, label, col in GROUPS:
    g = (dsid, reg)
    rs = recs.get(g, [])
    st = status.get(g, {})
    if not rs:
        print(f"{dsid+' / '+reg:<34}{0:>5}{st.get('skipped',0):>6}{st.get('error',0):>5}   (no results)")
        continue
    n = len(rs)
    k = sum(r["sequences"] for r in rs)         # windowed test, p_session < 0.05
    frac_rhy = np.mean([r["frac_rhythmic"] for r in rs])
    osc = np.array([r.get("oscillation_score", np.nan) for r in rs])   # frac sig windows
    # is the oscillation score systematically above the 0.05 chance level?
    _, p = stats.wilcoxon(osc - 0.05) if n >= 6 and np.isfinite(osc).all() else (np.nan, np.nan)
    print(f"{dsid+' / '+reg:<34}{n:>5}{st.get('skipped',0):>6}{st.get('error',0):>5}   "
          f"{k}/{n} = {100*k/n:>3.0f}%   "
          f"median osc-score={np.nanmedian(osc):.2f}  (Wilcoxon vs 0.05 p={p:.3g})"
          f"   mean rhythmic={100*frac_rhy:.0f}%")
    rows.append(dict(label=label, color=col, n=n, k=k, frac=k / n,
                     osc=osc, frac_rhy=frac_rhy, p=p))

if not rows:
    raise SystemExit("no results yet")

# ---------------- figure ----------------
fig, axes = plt.subplots(1, 2, figsize=(16, 6.6))

ax = axes[0]
x = np.arange(len(rows))
# hatch the groups whose 'sequences' are a known behavioural/task confound:
# V1 (lap cycle) and macaque EC (task-engagement blocks). See 05/10 controls.
CONFOUNDED = {"DANDI 001701\nvisual cortex\n(region control)",
              "DANDI 000897\nMACAQUE EC\nmental nav."}
bars = ax.bar(x, [100 * r["frac"] for r in rows], color=[r["color"] for r in rows],
              edgecolor="k", lw=0.6)
for b, r in zip(bars, rows):
    if r["label"] in CONFOUNDED:
        b.set_hatch("//")
ax.set_xticks(x)
ax.set_xticklabels([r["label"] for r in rows], fontsize=7.6, rotation=30, ha="right")
ax.set_ylabel("% of sessions with periodic sequences")
ax.set_ylim(0, 100)
ax.axhline(100 * 15 / 27, color="crimson", ls="--", lw=1.2,
           label="paper: 15/27 MEC wheel sessions (56%)")
from matplotlib.patches import Patch
confound_proxy = Patch(facecolor="white", edgecolor="k", hatch="///",
                       label="behavioural/task confound (see controls)")
h, _ = ax.get_legend_handles_labels()
ax.legend(handles=h + [confound_proxy], fontsize=8.5)
for xi, r in zip(x, rows):
    ax.text(xi, 100 * r["frac"] + 2, f"{r['k']}/{r['n']}", ha="center", fontsize=8.5)
ax.set_title("(a) Fraction of SESSIONS with periodic sequences\n"
             "(the quantity the paper reports — robust to session variability)",
             fontsize=10.5, pad=10)

ax = axes[1]
parts = ax.violinplot([r["osc"][np.isfinite(r["osc"])] for r in rows],
                      positions=x, showmedians=True, widths=0.8)
for pc, r in zip(parts["bodies"], rows):
    pc.set_facecolor(r["color"])
    pc.set_alpha(0.75)
ax.axhline(0.05, color="crimson", ls="--", lw=1.2, label="chance level (0.05)")
ax.set_xticks(x)
ax.set_xticklabels([r["label"] for r in rows], fontsize=7.6, rotation=30, ha="right")
ax.set_ylabel("oscillation score (fraction of significant windows)")
ax.legend(fontsize=8.5)
ax.set_title("(b) Per-session oscillation score\n(fraction of 300 s windows with sequences)",
             fontsize=10.5, pad=10)

fig.suptitle("Every session, every dandiset — minute-scale oscillatory sequences",
             fontsize=13.5, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.93])
out = f"{ROOT}/figures/SUMMARY_all_sessions.png"
fig.savefig(out, dpi=125)
print(f"\n-> {out}")
