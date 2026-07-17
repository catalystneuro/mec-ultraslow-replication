"""Is the CA1 null just insufficient power?

The honest worry about any null: the CA1 rest sessions (000552) have ~1.9x fewer
units than the wheel positive control (215 vs 404) and ~2.5x lower firing rates
(1.25 vs 3.18 Hz median). Either could hide a real effect. A null is only worth
reporting if the detector still fires under the SAME handicap.

So: degrade the positive control -- the paper's own wheel-in-darkness data, where
the effect is known present -- down to CA1's conditions, and ask whether the
sequences survive.

  (1) UNIT SUBSAMPLING: random subsets of n units, n = 50 .. all.
  (2) RATE MATCHING: independently thin each spike train to CA1's median rate.
      Thinning a rate-modulated process preserves the modulation and adds Poisson
      noise, so it reproduces CA1's SNR without inventing a different signal.
  (3) BOTH: n = 215 AND thinned to 1.25 Hz -- the fully CA1-matched control.

If the wheel still shows sequences under (3), then unit count and firing rate do
not explain the CA1 null. If it does not, the CA1 null is uninformative and must
be withdrawn.

usage: python3 scripts/13_power_analysis.py
"""

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
N_SHIFT = 200
N_REPEATS = 5

# measured from 000552 (see 11_ca1_rest_000552.py)
CA1_N_UNITS = 215
CA1_MEDIAN_RATE = 1.25
CA1_N_RANGE = (102, 221)


def thin(spike_times, p, rng):
    """Keep each spike independently with probability p."""
    return [s[rng.random(len(s)) < p] for s in spike_times]


def run(spike_times, t0, t1, seed):
    z = mo.build_activity(spike_times, t0, t1)[1][::DECIM]
    win = mo.windowed_sequence_test(z, mo.BIN_SIZE * DECIM,
                                    n_surrogates=N_SHIFT, seed=seed)
    return win["sequences"], win["frac_sig_windows"], win["p_session"]


def main():
    out = {}
    for mouse in ("104638", "102335"):
        s = loaders.load_ebrains(mouse, "Wheel-HeadFixed")
        rates = np.array([len(x) / s.duration for x in s.spike_times])
        med = float(np.median(rates))
        print(f"\n{'='*78}\n{mouse} wheel: {s.n_units} units, {s.duration/60:.1f} min, "
              f"median {med:.2f} Hz\n{'='*78}", flush=True)

        # ---- (1) unit subsampling, full rate
        grid = [n for n in (50, 75, 100, 150, 215, 275, 350, 404) if n <= s.n_units]
        curve = []
        for n in grid:
            hits, scores = [], []
            for r in range(N_REPEATS):
                rng = np.random.default_rng(1000 + r)
                idx = rng.choice(s.n_units, n, replace=False)
                seq, sc, _ = run([s.spike_times[i] for i in idx], s.t0, s.t1, seed=r)
                hits.append(seq); scores.append(sc)
            curve.append(dict(n=n, detect=float(np.mean(hits)),
                              osc=float(np.median(scores))))
            print(f"  n={n:>3}  detected {sum(hits)}/{N_REPEATS}  "
                  f"median osc-score={np.median(scores):.2f}", flush=True)

        # ---- (2)+(3) CA1-matched: n=215 AND thinned to CA1's rate
        p_keep = min(1.0, CA1_MEDIAN_RATE / med)
        matched = []
        n_match = min(CA1_N_UNITS, s.n_units)
        for r in range(N_REPEATS):
            rng = np.random.default_rng(2000 + r)
            idx = rng.choice(s.n_units, n_match, replace=False)
            st = thin([s.spike_times[i] for i in idx], p_keep, rng)
            seq, sc, p = run(st, s.t0, s.t1, seed=r)
            matched.append(dict(seq=bool(seq), osc=float(sc), p=float(p)))
        det = sum(m["seq"] for m in matched)
        print(f"  CA1-MATCHED (n={n_match}, thinned to {CA1_MEDIAN_RATE} Hz, "
              f"keep p={p_keep:.2f}): detected {det}/{N_REPEATS}, "
              f"median osc-score={np.median([m['osc'] for m in matched]):.2f}", flush=True)

        out[mouse] = dict(n_units=s.n_units, median_rate=med, curve=curve,
                          matched=matched, p_keep=p_keep, n_match=n_match)

    with open(f"{ROOT}/results/power_analysis.json", "w") as fh:
        json.dump(out, fh, indent=2)

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
    ax = axes[0]
    for mouse, col in (("104638", "#2166ac"), ("102335", "#92c5de")):
        d = out[mouse]
        ax.plot([c["n"] for c in d["curve"]], [100 * c["detect"] for c in d["curve"]],
                "o-", color=col, label=f"wheel {mouse} ({d['n_units']} units)")
    ax.axvspan(*CA1_N_RANGE, color="#5aae61", alpha=0.18,
               label=f"CA1 rest unit range ({CA1_N_RANGE[0]}-{CA1_N_RANGE[1]})")
    ax.axhline(0, color="grey", lw=0.6, ls=":")
    ax.set_xlabel("number of simultaneously recorded units (random subsets)")
    ax.set_ylabel("% of subsets where sequences are detected")
    ax.set_ylim(-6, 106)
    ax.legend(fontsize=8.5)
    ax.set_title("(a) Detection vs unit count, positive control\n"
                 "does the effect survive down to CA1's n?", fontsize=10.5)

    ax = axes[1]
    labels, vals, cols = [], [], []
    for mouse, col in (("104638", "#2166ac"), ("102335", "#92c5de")):
        d = out[mouse]
        labels.append(f"wheel {mouse}\nCA1-matched\n(n={d['n_match']}, {CA1_MEDIAN_RATE} Hz)")
        vals.append(100 * np.mean([m["seq"] for m in d["matched"]]))
        cols.append(col)
    labels.append(f"CA1 rest\n(actual)")
    vals.append(0.0)
    cols.append("#5aae61")
    ax.bar(range(len(vals)), vals, color=cols, edgecolor="k", lw=0.6)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("% with sequences detected")
    ax.set_ylim(0, 106)
    ax.set_title("(b) Positive control degraded to CA1's units AND firing rate\n"
                 "vs CA1's actual result", fontsize=10.5)

    fig.suptitle("Power analysis: is the CA1 null explained by too few / too sparse cells?",
                 fontsize=13, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    p = f"{ROOT}/figures/power_analysis.png"
    fig.savefig(p, dpi=125)
    print(f"\n-> {p}")


if __name__ == "__main__":
    main()
