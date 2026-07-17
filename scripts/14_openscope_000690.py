"""MEC vs its own negative-control regions, recorded simultaneously (DANDI 000690).

Allen OpenScope Vision2Hippocampus was never screened by this repo -- its title says
"Vision2Hippocampus" and its metadata never mentions MEC -- yet it holds the best
medial entorhinal coverage on DANDI: up to 726 good MEC units in one session,
against 404 and 273 in the paper's own Neuropixels recordings.

What it uniquely offers: the paper's Fig. 5 design WITHIN a session. The original
compared MEC (15/27 sessions) against parasubiculum (0/25) and visual cortex (0/19)
across different mice, days and brain states. Here all three sit on the same probes
in the same animal at the same moment, under identical input.

That matters because of this panel's recurring failure mode: apparent sequences that
turn out to be behavioural or task structure (V1's lap cycle in 001701, task blocks
in the macaque). Here the stimulus is common to all regions by construction, so a
MEC-vs-V1 difference cannot be manufactured by stimulus structure -- if it could, V1
would show it at least as strongly, being the visual area.

UNIT-COUNT MATCHING. MEC yields ~2.4x more units than V1 in the same session (668 vs
276), so an unmatched contrast would confound region with statistical power. Each
session is therefore subsampled to a common n across its regions (capped at
MAX_UNITS). 13_power_analysis.py licenses this: on the positive control the effect is
detected 5/5 down to n=50, so matching at ~250 costs nothing detectable.

What this CANNOT do: test the paper's condition. There is no darkness and no
task-free block. A null is ambiguous between "sequences reset under sensory drive"
(the authors' own speculation) and "no effect"; only the relative region contrast is
robust.

usage: python3 scripts/14_openscope_000690.py [--limit N]
"""

import argparse
import json
import os
import re
import sys
import time
import traceback

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import loaders
import mec_oscillation as mo

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f"{ROOT}/results/sessions"
os.makedirs(OUT, exist_ok=True)

N_ISI = 50
N_SHIFT = 200
DECIM = 8
MIN_UNITS = 20
MIN_MINUTES = 10
MAX_UNITS = 250          # cap; see UNIT-COUNT MATCHING above
REGIONS = ("MEC", "V1", "PaS")


def key(asset, region):
    safe = re.sub(r"[^A-Za-z0-9]+", "-", asset).strip("-")
    return f"000690__{safe}__{region}"


def analyze(spike_times, t0, t1):
    zt = mo.build_activity(spike_times, t0, t1)[1]
    zd = zt[::DECIM]
    bin_s = mo.BIN_SIZE * DECIM
    sc = mo.single_cell_test(spike_times, t0, t1, n_surrogates=N_ISI)
    pop = mo.population_sequence_test(zd, n_surrogates=N_SHIFT)
    win = mo.windowed_sequence_test(zd, bin_s, n_surrogates=N_SHIFT)
    sig = sc["sig"]
    return dict(
        frac_rhythmic=float(sig.mean()), n_rhythmic=int(sig.sum()),
        median_peak_hz=float(np.median(sc["f_peak"][sig])) if sig.any() else None,
        pc12_var=pop["pc12_var"], pc12_z=pop["pc12_z"], p_pc12=pop["p_pc12"],
        rot_z=pop["rot_z"], p_rotation=pop["p_rotation"],
        n_windows=win["n_windows"], n_sig_windows=win["n_sig_windows"],
        oscillation_score=win["frac_sig_windows"], p_session=win["p_session"],
        best_window_rot_z=win["best_rot_z"], sequences=bool(win["sequences"]),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--retry", action="store_true")
    args = ap.parse_args()

    paths = loaders.OPENSCOPE_MEC_SESSIONS[: args.limit] if args.limit \
        else loaders.OPENSCOPE_MEC_SESSIONS
    print(f"{'#' * 90}\n# 000690: {len(paths)} sessions x {len(REGIONS)} regions "
          f"(unit-matched within session)\n{'#' * 90}", flush=True)

    for i, asset in enumerate(paths):
        # load every region first so they can be matched on unit count
        sess = {}
        for reg in REGIONS:
            try:
                s = loaders.load_000690(asset, region_select=reg)
                if s.n_units >= MIN_UNITS and s.duration / 60 >= MIN_MINUTES:
                    sess[reg] = s
            except Exception:
                pass
        if not sess:
            print(f"[{i+1}/{len(paths)}] {asset[:40]}  no usable regions", flush=True)
            continue

        n_common = min(min(s.n_units for s in sess.values()), MAX_UNITS)
        rng = np.random.default_rng(0)

        for reg, s in sess.items():
            k = key(asset, reg)
            p = f"{OUT}/{k}.json"
            if os.path.exists(p) and not args.retry:
                print(f"[{i+1}/{len(paths)}] {k[:66]}  (cached)", flush=True)
                continue
            t0 = time.time()
            rec = dict(dandiset="000690", asset=asset, region_select=reg)
            try:
                idx = rng.choice(s.n_units, n_common, replace=False)
                st = [s.spike_times[j] for j in idx]
                rec.update(status="ok", n_units=int(n_common),
                           n_units_available=int(s.n_units),
                           duration_min=float(s.duration / 60), species=s.species,
                           region=s.region, description=s.description,
                           subject=s.subject)
                rec.update(analyze(st, s.t0, s.t1))
            except Exception as e:
                rec.update(status="error", error=f"{type(e).__name__}: {e}",
                           traceback=traceback.format_exc()[-1200:])
            with open(p, "w") as fh:
                json.dump(rec, fh, indent=2)
            dt = time.time() - t0
            if rec["status"] == "ok":
                msg = (f"{reg:3s} units={rec['n_units']:>3}/{rec['n_units_available']:<3} "
                       f"{rec['duration_min']:>5.1f}min rhythmic={100*rec['frac_rhythmic']:>3.0f}% "
                       f"rot_z={rec['rot_z']:+5.1f} seq={'YES' if rec['sequences'] else 'no '}")
            else:
                msg = f"{reg:3s} ERROR {rec['error'][:50]}"
            print(f"[{i+1}/{len(paths)}] {msg}  [{dt:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
