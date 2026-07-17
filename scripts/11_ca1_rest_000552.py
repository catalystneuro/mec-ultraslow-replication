"""Ask the paper's question of hippocampal CA1 during long TASK-FREE rest
(dandiset 000552, Huszar et al.).

Why: the panel so far tests MEC only under navigation, because every entorhinal
dataset on DANDI runs a task from start to finish. 000552 inverts that trade --
it gives up the region (CA1, not MEC) to buy the condition (no task, no rewards,
minimal sensory drive, 39-148 min continuous). It is the only public data that
can be asked the question in a task-free state.

What it can and cannot show. In the paper's own logic this is a REGION control,
so a null is the EXPECTED result and corroborates rather than tests: MEC carries
a low-dimensional code with ~94% of cells locked to the sequence, whereas
hippocampus carries a high-dimensional code where reported sequences involve ~5%
of the network. What the null does add is a condition the authors flag as open
("whether the ultraslow oscillatory sequences are present across a broader
spectrum of behaviours, including sleep and free exploration").

Selection (see scan in the docstring of loaders.load_000552): sessions with
>=100 units and >=30 min of continuous task-free Awake rest. 24 of 64 qualify.

usage:
  python3 scripts/11_ca1_rest_000552.py            # top sessions, resuming
  python3 scripts/11_ca1_rest_000552.py --limit 3
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

# Sessions with >=100 units and >=30 min continuous task-free Awake rest,
# ranked by units x duration. Determined by streaming every session's epochs,
# trials and SleepStates tables (24 of 64 processed sessions qualify).
SESSIONS = [
    "sub-e15-13f1/sub-e15-13f1_ses-e15-13f1-220119_behavior+ecephys.nwb",
    "sub-e14-2m2/sub-e14-2m2_ses-e14-2m2-201013_behavior+ecephys.nwb",
    "sub-e15-13f1/sub-e15-13f1_ses-e15-13f1-220117_behavior+ecephys.nwb",
    "sub-e15-13f1/sub-e15-13f1_ses-e15-13f1-220118_behavior+ecephys.nwb",
    "sub-e16-3m1/sub-e16-3m1_ses-e16-3m1-210111_behavior+ecephys.nwb",
    "sub-e16-3m1/sub-e16-3m1_ses-e16-3m1-210201_behavior+ecephys.nwb",
]


def key(asset, region):
    safe = re.sub(r"[^A-Za-z0-9]+", "-", asset).strip("-")
    return f"000552__{safe}__{region}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--retry", action="store_true")
    ap.add_argument("--state", default="Awake")
    args = ap.parse_args()

    paths = SESSIONS[: args.limit] if args.limit else SESSIONS
    region = f"CA1-{args.state}"
    print(f"{'#' * 90}\n# 000552: {len(paths)} task-free {args.state} rest sessions\n{'#' * 90}",
          flush=True)

    for i, asset in enumerate(paths):
        k = key(asset, region)
        p = f"{OUT}/{k}.json"
        if os.path.exists(p) and not args.retry:
            print(f"[{i+1}/{len(paths)}] {k[:72]}  (cached)", flush=True)
            continue

        t0 = time.time()
        rec = dict(dandiset="000552", asset=asset, region_select=region)
        try:
            sess = loaders.load_000552(asset, state=args.state)
            if sess.n_units < MIN_UNITS or sess.duration / 60 < MIN_MINUTES:
                rec.update(status="skipped", n_units=int(sess.n_units),
                           duration_min=float(sess.duration / 60),
                           reason=f"n_units<{MIN_UNITS} or dur<{MIN_MINUTES}min")
            else:
                zt = mo.build_activity(sess.spike_times, sess.t0, sess.t1)[1]
                zd = zt[::DECIM]
                bin_s = mo.BIN_SIZE * DECIM
                sc = mo.single_cell_test(sess.spike_times, sess.t0, sess.t1,
                                         n_surrogates=N_ISI)
                pop = mo.population_sequence_test(zd, n_surrogates=N_SHIFT)
                win = mo.windowed_sequence_test(zd, bin_s, n_surrogates=N_SHIFT)
                sig = sc["sig"]
                rec.update(
                    status="ok",
                    n_units=int(sess.n_units),
                    duration_min=float(sess.duration / 60),
                    species=sess.species, region=sess.region,
                    description=sess.description, subject=sess.subject,
                    frac_rhythmic=float(sig.mean()), n_rhythmic=int(sig.sum()),
                    median_peak_hz=(float(np.median(sc["f_peak"][sig]))
                                    if sig.any() else None),
                    pc12_var=pop["pc12_var"], pc12_z=pop["pc12_z"],
                    p_pc12=pop["p_pc12"], rot_z=pop["rot_z"],
                    p_rotation=pop["p_rotation"],
                    n_windows=win["n_windows"], n_sig_windows=win["n_sig_windows"],
                    oscillation_score=win["frac_sig_windows"],
                    p_session=win["p_session"],
                    best_window_rot_z=win["best_rot_z"],
                    sequences=bool(win["sequences"]),
                )
        except Exception as e:
            rec.update(status="error", error=f"{type(e).__name__}: {e}",
                       traceback=traceback.format_exc()[-1500:])

        with open(p, "w") as fh:
            json.dump(rec, fh, indent=2)

        dt = time.time() - t0
        if rec["status"] == "ok":
            msg = (f"units={rec['n_units']:>3} {rec['duration_min']:>5.1f}min "
                   f"rhythmic={100*rec['frac_rhythmic']:>3.0f}% "
                   f"rot_z={rec['rot_z']:+5.1f} "
                   f"seq={'YES' if rec['sequences'] else 'no '}")
        elif rec["status"] == "skipped":
            msg = f"SKIP ({rec['reason']})"
        else:
            msg = f"ERROR {rec['error'][:60]}"
        print(f"[{i+1}/{len(paths)}] {k[:56]:<56} {msg}  [{dt:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
