"""Is MEC's ultraslow oscillation just arousal? The control the paper cannot run.

The paper names this alternative itself: the "widespread nature of the ultraslow
oscillatory activity in individual neurons would be consistent with a role for
ascending neuromodulatory AROUSAL-associated brain-stem circuits in controlling
these oscillations", citing infraslow PUPIL oscillation (Blasiak 2013) and vascular
contributions to ultra-slow signals (Drew 2020). And 0.006-0.008 Hz sits squarely in
the arousal/vasomotion range.

But the paper has no pupil, no vascular and no independent arousal measure anywhere.
The control is orphaned.

DANDI 000690 can run it: EyeTracking (pupil area) and running speed alongside the
archive's best MEC coverage. Sequences cannot be tested there -- the screen never goes
dark -- but sequences are not what the arousal hypothesis is about. The SINGLE-CELL
oscillation is, and that is measurable here.

Method follows 05_behavior_control.py: run the SAME detector (spectrum of the
z-scored autocorrelogram) on the arousal signal as on the neural population, and ask
whether their ultraslow peaks coincide. Coupling is then tested with band coherence
against an independent CIRCULAR SHIFT of the arousal trace -- the same null as
population_sequence_test, which preserves each signal's own spectrum and destroys
only their relative timing.

That null is essential and a naive statistic fails here. Both signals are 1/f-like,
so their LOG SPECTRA correlate ~0.94 whether or not they share an oscillation; on the
first two sessions that statistic reported near-perfect agreement and is meaningless.
It is the same trap as "band power cannot detect the effect" in the README.

Reading the result:
  - coherence above the shift null, peaks coinciding -> the ultraslow oscillation
    tracks arousal, supporting the authors' own ascending-arousal account and leaving
    MEC-specific SEQUENCE formation as the residual claim.
  - coherence at the null -> the oscillation is more intrinsic than they guessed.
Neither outcome bears on the sequence claim; this is about the oscillation.

usage: python3 scripts/15_arousal_control.py [--limit N]
"""

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal as sps
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(__file__))
import loaders
import mec_oscillation as mo

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DECIM = 8
BAND = mo.BAND
N_SHIFT = 200


def band_coherence(a, b, fs, nperseg=1024):
    """Mean magnitude-squared coherence in the ultraslow band.

    Coherence is biased upward with few segments, but the circular-shift null is
    computed identically and carries the same bias, so the comparison is valid.
    """
    n = min(len(a), len(b), nperseg)
    f, c = sps.coherence(a, b, fs=fs, nperseg=n, noverlap=n // 2)
    m = (f >= BAND[0]) & (f <= BAND[1])
    return float(np.nanmean(c[m])) if m.any() else np.nan


def resample_to(t_src, v_src, t_dst):
    ok = np.isfinite(t_src) & np.isfinite(v_src)
    if ok.sum() < 10:
        return None
    return np.interp(t_dst, t_src[ok], v_src[ok])


def spectrum_of(x, fs):
    """Identical detector to the neural data (05_behavior_control.py)."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 100 or x.std() == 0:
        return None, None
    return mo.acg_spectrum(x - x.mean(), fs=fs)


def band_peak(f, p):
    m = (f >= BAND[0]) & (f <= BAND[1])
    if not m.any():
        return np.nan, np.nan
    i = np.argmax(p[m])
    return f[m][i], p[m][i]


def analyse(asset):
    s = loaders.load_000690(asset, region_select="MEC")
    z = mo.build_activity(s.spike_times, s.t0, s.t1)[1]
    bin_s = mo.BIN_SIZE * DECIM
    zd = z[::DECIM]
    t = s.t0 + np.arange(zd.shape[0]) * bin_s
    fs = 1.0 / bin_s

    # population signal: mean smoothed activity (the "global" ultraslow drive)
    pop = zd.mean(axis=1)
    pop = gaussian_filter1d(pop, 1.0)

    tp, pupil, tr, run = loaders.load_000690_arousal(asset)
    out = dict(asset=asset, subject=s.subject, n_units=s.n_units,
               duration_min=float(s.duration / 60))

    fpop, ppop = spectrum_of(pop, fs)
    out["pop_peak_hz"] = float(band_peak(fpop, ppop)[0]) if fpop is not None else None

    for name, tt, vv in (("pupil", tp, pupil), ("running", tr, run)):
        if vv is None:
            out[f"{name}_peak_hz"] = None
            continue
        r = resample_to(np.asarray(tt, float), np.asarray(vv, float), t)
        if r is None:
            out[f"{name}_peak_hz"] = None
            continue
        r = gaussian_filter1d(r, 5.0 / bin_s)
        f_, p_ = spectrum_of(r, fs)
        out[f"{name}_peak_hz"] = float(band_peak(f_, p_)[0]) if f_ is not None else None

        # NOTE: a correlation between the two LOG SPECTRA is not usable here. Both
        # signals are 1/f-like, so any two of them correlate ~0.9 in log space
        # whether or not they share an oscillation -- the same trap documented in
        # "Two statistical traps found and fixed". Coupling must be tested against a
        # null that preserves each signal's own spectrum and destroys only their
        # relative timing: an independent CIRCULAR SHIFT of the arousal trace, the
        # same null used by population_sequence_test.
        ok = np.isfinite(r) & np.isfinite(pop)
        obs_r = float(np.corrcoef(pop[ok], r[ok])[0, 1])
        obs_c = float(band_coherence(pop[ok], r[ok], fs))
        rng = np.random.default_rng(0)
        T = ok.sum()
        null_r, null_c = [], []
        for _ in range(N_SHIFT):
            sh = np.roll(r[ok], int(rng.integers(T)))
            null_r.append(abs(np.corrcoef(pop[ok], sh)[0, 1]))
            null_c.append(band_coherence(pop[ok], sh, fs))
        null_r, null_c = np.array(null_r), np.array(null_c)
        out[f"{name}_trace_corr"] = obs_r
        out[f"{name}_trace_p"] = float((1 + (null_r >= abs(obs_r)).sum()) / (1 + N_SHIFT))
        out[f"{name}_coherence"] = obs_c
        out[f"{name}_coh_null"] = float(null_c.mean())
        out[f"{name}_coh_p"] = float((1 + (null_c >= obs_c).sum()) / (1 + N_SHIFT))
        out[f"{name}_coh_z"] = float((obs_c - null_c.mean()) / (null_c.std() + 1e-12))
        out[f"_{name}_trace"] = r
    out["_pop"] = pop
    out["_t"] = t
    out["_fpop"] = fpop
    out["_ppop"] = ppop
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    paths = loaders.OPENSCOPE_MEC_SESSIONS[: args.limit] if args.limit \
        else loaders.OPENSCOPE_MEC_SESSIONS

    rows = []
    for i, a in enumerate(paths):
        try:
            r = analyse(a)
            rows.append(r)
            print(f"[{i+1}/{len(paths)}] {r['subject']:<10} units={r['n_units']:>3} "
                  f"pop_peak={r['pop_peak_hz']:.4f}Hz "
                  f"pupil_peak={r.get('pupil_peak_hz') and round(r['pupil_peak_hz'],4)} "
                  f"coh={r.get('pupil_coherence', float('nan')):.3f} "
                  f"null={r.get('pupil_coh_null', float('nan')):.3f} "
                  f"z={r.get('pupil_coh_z', float('nan')):+.1f} "
                  f"p={r.get('pupil_coh_p', float('nan')):.3f}", flush=True)
        except Exception as e:
            print(f"[{i+1}/{len(paths)}] {a[:40]} ERR {type(e).__name__}: {e}", flush=True)

    clean = [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]
    json.dump(clean, open(f"{ROOT}/results/arousal_000690.json", "w"), indent=2)

    # ---------------- figure ----------------
    ex = rows[0]
    fig, axes = plt.subplots(2, 2, figsize=(15, 8.5))

    ax = axes[0, 0]
    tt = (ex["_t"] - ex["_t"][0]) / 60
    ax.plot(tt, (ex["_pop"] - ex["_pop"].mean()) / ex["_pop"].std(), lw=0.7,
            label="MEC population (z)", color="#2166ac")
    if "_pupil_trace" in ex:
        p = ex["_pupil_trace"]
        ax.plot(tt, (p - np.nanmean(p)) / np.nanstd(p), lw=0.7, alpha=0.8,
                label="pupil area (z)", color="#d6604d")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("z")
    ax.legend(fontsize=8)
    ax.set_title(f"(a) {ex['subject']}: MEC population vs pupil", fontsize=10.5)

    ax = axes[0, 1]
    if ex["_fpop"] is not None:
        m = (ex["_fpop"] >= BAND[0]) & (ex["_fpop"] <= BAND[1])
        ax.semilogy(ex["_fpop"][m], ex["_ppop"][m], color="#2166ac", label="MEC population")
    if "_pupil_trace" in ex:
        f_, p_ = spectrum_of(ex["_pupil_trace"], 1.0 / (mo.BIN_SIZE * DECIM))
        if f_ is not None:
            m = (f_ >= BAND[0]) & (f_ <= BAND[1])
            ax.semilogy(f_[m], p_[m], color="#d6604d", label="pupil")
    ax.axvspan(0.006, 0.008, color="green", alpha=0.2, label="paper's 0.006-0.008 Hz")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("PSD of ACG (a.u.)")
    ax.legend(fontsize=8)
    ax.set_title("(b) Same detector on both signals", fontsize=10.5)

    ax = axes[1, 0]
    pk_pop = [r["pop_peak_hz"] for r in rows if r.get("pupil_peak_hz")]
    pk_pup = [r["pupil_peak_hz"] for r in rows if r.get("pupil_peak_hz")]
    ax.scatter(pk_pup, pk_pop, s=80, color="#2166ac", edgecolor="k")
    lim = [0, max(pk_pop + pk_pup + [0.02]) * 1.1]
    ax.plot(lim, lim, "k--", lw=1, label="identity (arousal = neural)")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("pupil ultraslow peak (Hz)")
    ax.set_ylabel("MEC population peak (Hz)")
    ax.legend(fontsize=8)
    ax.set_title("(c) Do the peaks coincide? (each point = session)", fontsize=10.5)

    ax = axes[1, 1]
    tc = [abs(r["pupil_trace_corr"]) for r in rows if "pupil_trace_corr" in r]
    co = [r["pupil_coherence"] for r in rows if "pupil_coherence" in r]
    cn = [r["pupil_coh_null"] for r in rows if "pupil_coh_null" in r]
    ax.boxplot([co, cn, tc], tick_labels=["coherence\nobserved", "coherence\nshift null",
                                          "|trace corr|"])
    ax.axhline(0, color="grey", lw=0.6, ls=":")
    ax.set_ylabel("MEC-pupil coupling")
    ax.set_title(f"(d) MEC-pupil coupling across {len(tc)} sessions", fontsize=10.5)

    fig.suptitle("Is MEC's ultraslow oscillation arousal? (DANDI 000690 — the control "
                 "the original paper could not run)", fontsize=13, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    p = f"{ROOT}/figures/arousal_control_000690.png"
    fig.savefig(p, dpi=125)
    print(f"\n-> {p}")

    if tc:
        coh = [r["pupil_coherence"] for r in rows if "pupil_coherence" in r]
        cohn = [r["pupil_coh_null"] for r in rows if "pupil_coh_null" in r]
        cohp = [r["pupil_coh_p"] for r in rows if "pupil_coh_p" in r]
        cohz = [r["pupil_coh_z"] for r in rows if "pupil_coh_z" in r]
        d = [abs(a - b) for a, b in zip(pk_pop, pk_pup)]
        print(f"\nSUMMARY across {len(tc)} sessions (MEC population vs PUPIL):")
        print(f"  |peak freq difference| median   = {np.median(d):.4f} Hz")
        print(f"  band coherence  observed median = {np.median(coh):.3f}")
        print(f"  band coherence  shift-null mean = {np.median(cohn):.3f}")
        print(f"  coherence z vs null      median = {np.median(cohz):+.1f}")
        print(f"  sessions with coherence p<0.05  = {sum(1 for x in cohp if x < 0.05)}/{len(cohp)}")
        print(f"  |trace corr| median             = {np.median(tc):.2f}")


if __name__ == "__main__":
    main()
