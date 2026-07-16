"""Run the validated pipeline over EVERY session in each dandiset.

Motivation: the paper found oscillatory sequences in only 15 of 27 MEC wheel
sessions (Fig. 5g) -- the effect is strongly session-variable even under ideal
conditions. A single-session null is therefore weak evidence. The quantity the
paper actually reports is the FRACTION OF SESSIONS with sequences, so that is
what we compute.

RESUMABLE. One JSON per session under results/sessions/. A session that already
has a result file is skipped, so the run can be stopped and restarted freely.
Failures are recorded too (status='error') and are retried only with --retry.

usage:
  python3 07_run_all_sessions.py                  # everything, resuming
  python3 07_run_all_sessions.py --only 001701    # one dandiset
  python3 07_run_all_sessions.py --retry          # re-attempt failed sessions
  python3 07_run_all_sessions.py --limit 5        # first N sessions per dandiset
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

N_ISI = 50        # single-cell test (non-specific; kept for reporting only)
N_SHIFT = 200     # population sequence test -- the load-bearing one
DECIM = 8
MIN_UNITS = 20    # need a population to speak of a population sequence
MIN_MINUTES = 10  # need several cycles of a minute-scale rhythm


def key(dsid, asset, region):
    safe = re.sub(r"[^A-Za-z0-9]+", "-", f"{asset}").strip("-")
    return f"{dsid}__{safe}__{region}"


def already_done(k, retry):
    p = f"{OUT}/{k}.json"
    if not os.path.exists(p):
        return False
    if retry:
        try:
            if json.load(open(p)).get("status") == "error":
                return False
        except Exception:
            return False
    return True


def analyze_session(sess):
    """Run all tests on one loaded session. Returns the result dict.

    'sequences' uses the WINDOWED test, because the whole-session rotation index
    washes out intermittent sequences (it gave a false negative on the paper's
    own mouse 102335). The detector is calibrated to the strongest positive
    control, mouse 104638 (windowed p ~ 1e-6). Mouse 102335 -- the harder of the
    paper's two Neuropixels sessions, whose sequences the paper itself calls
    subtle and subset-dependent -- does NOT reach significance here, and we
    deliberately do not tune to force it (that would inflate false positives). So
    a null below means "no 104638-strength sequences"; weaker intermittent
    structure can be missed.
    """
    zt = mo.build_activity(sess.spike_times, sess.t0, sess.t1)[1]
    zd = zt[::DECIM]
    bin_s = mo.BIN_SIZE * DECIM

    sc = mo.single_cell_test(sess.spike_times, sess.t0, sess.t1, n_surrogates=N_ISI)
    pop = mo.population_sequence_test(zd, n_surrogates=N_SHIFT)
    win = mo.windowed_sequence_test(zd, bin_s, n_surrogates=N_SHIFT)

    sig = sc["sig"]
    return dict(
        status="ok",
        n_units=int(sess.n_units),
        duration_min=float(sess.duration / 60),
        species=sess.species,
        region=sess.region,
        description=sess.description,
        subject=sess.subject,
        frac_rhythmic=float(sig.mean()),
        n_rhythmic=int(sig.sum()),
        median_peak_hz=float(np.median(sc["f_peak"][sig])) if sig.any() else None,
        # whole-session population stats (kept for reference)
        pc12_var=pop["pc12_var"], pc12_z=pop["pc12_z"], p_pc12=pop["p_pc12"],
        rot_z=pop["rot_z"], p_rotation=pop["p_rotation"],
        # windowed sequence test (the primary criterion)
        n_windows=win["n_windows"], n_sig_windows=win["n_sig_windows"],
        oscillation_score=win["frac_sig_windows"], p_session=win["p_session"],
        best_window_rot_z=win["best_rot_z"],
        sequences=bool(win["sequences"]),
    )


def jobs_for(dsid, limit=None):
    """(asset_path, region, loader) for every analyzable session in a dandiset."""
    if dsid == "ebrains":
        out = []
        for mouse in ("104638", "102335"):
            for task in ("Wheel-HeadFixed", "OpenField"):
                out.append((f"{mouse}/{task}", task,
                            lambda m=mouse, tk=task: loaders.load_ebrains(m, tk)))
        return out

    if dsid == "000053":
        paths = [p for p in loaders.list_assets(dsid) if "ecephys" in p]
        regions = ["MEC"]
    elif dsid == "001701":
        paths = loaders.list_assets(dsid)
        regions = ["MEC", "V1"]     # V1 = within-session region control
    elif dsid == "000897":
        # each session appears twice: '..._behavior+ecephys.nwb' (has the units
        # table) and a bare '..._ecephys.nwb' (raw only, no units). Keep the former.
        paths = [p for p in loaders.list_assets(dsid) if "behavior+ecephys" in p]
        regions = ["EC"]
    else:
        raise ValueError(dsid)

    if limit:
        paths = paths[:limit]
    return [(p, r, (lambda pp=p, rr=r: loaders.load(dsid, pp, region_select=rr)))
            for p in paths for r in regions]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    ap.add_argument("--retry", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    dsids = [args.only] if args.only else ["ebrains", "000897", "001701", "000053"]

    for dsid in dsids:
        jobs = jobs_for(dsid, args.limit)
        print(f"\n{'#' * 90}\n# {dsid}: {len(jobs)} session x region jobs\n{'#' * 90}", flush=True)

        for i, (asset, region, loader) in enumerate(jobs):
            k = key(dsid, asset, region)
            if already_done(k, args.retry):
                print(f"[{i+1}/{len(jobs)}] {k[:72]}  (cached)", flush=True)
                continue

            t0 = time.time()
            rec = dict(dandiset=dsid, asset=asset, region_select=region)
            try:
                sess = loader()
                if sess.n_units < MIN_UNITS or sess.duration / 60 < MIN_MINUTES:
                    rec.update(status="skipped", n_units=int(sess.n_units),
                               duration_min=float(sess.duration / 60),
                               reason=f"n_units<{MIN_UNITS} or dur<{MIN_MINUTES}min")
                else:
                    rec.update(analyze_session(sess))
            except Exception as e:
                rec.update(status="error", error=f"{type(e).__name__}: {e}",
                           traceback=traceback.format_exc()[-1500:])

            with open(f"{OUT}/{k}.json", "w") as fh:
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
            print(f"[{i+1}/{len(jobs)}] {k[:60]:<60} {msg}  [{dt:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
