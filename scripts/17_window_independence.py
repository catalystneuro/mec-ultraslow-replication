"""Re-run the whole panel with overlapping AND independent windows, side by side.

The session-level criterion everywhere in this repo is
`binomtest(n_sig_windows, n_windows, 0.05)`. With the defaults (300 s window, 100 s
step) adjacent windows share 2/3 of their samples, so that binomial counts ~3x more
evidence than exists: it is anti-conservative. Non-overlapping windows are exactly
every 3rd window of the same grid, so both statistics come from one surrogate run
(`windowed_sequence_test` now returns both).

Why it matters, and in which direction:
  - POSITIVES weaken. On the paper's own mouse 102335 the overlapping test gives
    3/13, p=0.025 (significant); independent windows give 1/5, p=0.23 (not). That
    session's median observed rotation (0.498) is BELOW its median null (0.588).
  - NULLS strengthen, or are untouched. A liberal test that found nothing would find
    nothing under a stricter one.

So this re-run cannot rescue any null in the panel; it can only cost positives --
including one of the two positive controls.

Writes results/sessions_indep/*.json. Resumable. Windowed test only (the statistic
under review); single-cell and whole-session stats are unchanged and already stored.

usage: python3 scripts/17_window_independence.py [--only 000690] [--limit N]
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
OUT = f"{ROOT}/results/sessions_indep"
os.makedirs(OUT, exist_ok=True)

N_SHIFT = 200
DECIM = 8
MIN_UNITS = 20
MIN_MINUTES = 10
MAX_UNITS = 250


def key(dsid, asset, region):
    safe = re.sub(r"[^A-Za-z0-9]+", "-", f"{asset}").strip("-")
    return f"{dsid}__{safe}__{region}"


def score(spike_times, t0, t1):
    z = mo.build_activity(spike_times, t0, t1)[1][::DECIM]
    w = mo.windowed_sequence_test(z, mo.BIN_SIZE * DECIM, n_surrogates=N_SHIFT)
    return {k: w[k] for k in
            ("n_windows", "n_sig_windows", "frac_sig_windows", "p_session",
             "n_windows_indep", "n_sig_windows_indep", "frac_sig_windows_indep",
             "p_session_indep", "median_rot", "median_rot_null", "best_rot_z",
             "sequences", "sequences_indep")}


def jobs(dsid, limit=None):
    if dsid == "ebrains":
        return [(f"{m}/{t}", t, (lambda m=m, t=t: loaders.load_ebrains(m, t)))
                for m in ("104638", "102335") for t in ("Wheel-HeadFixed", "OpenField")]
    if dsid == "000552":
        return [(p, "CA1-Awake", (lambda p=p: loaders.load_000552(p)))
                for p in (["sub-e15-13f1/sub-e15-13f1_ses-e15-13f1-220119_behavior+ecephys.nwb",
                           "sub-e14-2m2/sub-e14-2m2_ses-e14-2m2-201013_behavior+ecephys.nwb",
                           "sub-e15-13f1/sub-e15-13f1_ses-e15-13f1-220117_behavior+ecephys.nwb",
                           "sub-e15-13f1/sub-e15-13f1_ses-e15-13f1-220118_behavior+ecephys.nwb",
                           "sub-e16-3m1/sub-e16-3m1_ses-e16-3m1-210111_behavior+ecephys.nwb",
                           "sub-e16-3m1/sub-e16-3m1_ses-e16-3m1-210201_behavior+ecephys.nwb"])[:limit]]
    if dsid == "000690":
        return [(p, r, None) for p in loaders.OPENSCOPE_MEC_SESSIONS[:limit]
                for r in ("MEC", "V1", "PaS")]
    if dsid == "000053":
        paths = [p for p in loaders.list_assets(dsid) if "ecephys" in p][:limit]
        regions = ["MEC"]
    elif dsid == "001701":
        paths = loaders.list_assets(dsid)[:limit]
        regions = ["MEC", "V1"]
    elif dsid == "000897":
        paths = [p for p in loaders.list_assets(dsid) if "behavior+ecephys" in p][:limit]
        regions = ["EC"]
    else:
        raise ValueError(dsid)
    return [(p, r, (lambda pp=p, rr=r: loaders.load(dsid, pp, region_select=rr)))
            for p in paths for r in regions]


def run_690(asset, limit=None):
    """000690 needs its regions unit-matched within session, as in 14_openscope."""
    sess = {}
    for r in ("MEC", "V1", "PaS"):
        try:
            s = loaders.load_000690(asset, region_select=r)
            if s.n_units >= MIN_UNITS and s.duration / 60 >= MIN_MINUTES:
                sess[r] = s
        except Exception:
            pass
    if not sess:
        return {}
    n = min(min(s.n_units for s in sess.values()), MAX_UNITS)
    rng = np.random.default_rng(0)
    out = {}
    for r, s in sess.items():
        idx = rng.choice(s.n_units, n, replace=False)
        out[r] = (s, [s.spike_times[j] for j in idx], n)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    dsids = [args.only] if args.only else \
        ["ebrains", "000552", "000690", "000053", "000897", "001701"]

    for dsid in dsids:
        J = jobs(dsid, args.limit)
        print(f"\n{'#'*88}\n# {dsid}: {len(J)} jobs\n{'#'*88}", flush=True)

        cache690 = {}
        for i, (asset, region, ld) in enumerate(J):
            k = key(dsid, asset, region)
            p = f"{OUT}/{k}.json"
            if os.path.exists(p):
                print(f"[{i+1}/{len(J)}] {k[:64]} (cached)", flush=True)
                continue
            t0 = time.time()
            rec = dict(dandiset=dsid, asset=asset, region_select=region)
            try:
                if dsid == "000690":
                    if asset not in cache690:
                        cache690 = {asset: run_690(asset)}
                    d = cache690[asset]
                    if region not in d:
                        rec.update(status="skipped", reason="region absent")
                    else:
                        s, st, n = d[region]
                        rec.update(status="ok", n_units=int(n), subject=s.subject,
                                   duration_min=float(s.duration / 60))
                        rec.update(score(st, s.t0, s.t1))
                else:
                    s = ld()
                    if s.n_units < MIN_UNITS or s.duration / 60 < MIN_MINUTES:
                        rec.update(status="skipped", n_units=int(s.n_units),
                                   duration_min=float(s.duration / 60), reason="too small")
                    else:
                        rec.update(status="ok", n_units=int(s.n_units), subject=s.subject,
                                   duration_min=float(s.duration / 60))
                        rec.update(score(s.spike_times, s.t0, s.t1))
            except Exception as e:
                rec.update(status="error", error=f"{type(e).__name__}: {e}",
                           traceback=traceback.format_exc()[-800:])
            json.dump(rec, open(p, "w"), indent=2)
            dt = time.time() - t0
            if rec["status"] == "ok":
                msg = (f"overlap {rec['n_sig_windows']:>2}/{rec['n_windows']:<3} "
                       f"p={rec['p_session']:.4f} {'SEQ' if rec['sequences'] else 'no '} | "
                       f"INDEP {rec['n_sig_windows_indep']:>2}/{rec['n_windows_indep']:<2} "
                       f"p={rec['p_session_indep']:.4f} "
                       f"{'SEQ' if rec['sequences_indep'] else 'no '}")
            else:
                msg = rec["status"] + " " + str(rec.get("reason", rec.get("error", "")))[:40]
            print(f"[{i+1}/{len(J)}] {k[:52]:<52} {msg} [{dt:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
