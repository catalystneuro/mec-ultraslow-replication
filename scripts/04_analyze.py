"""Definitive analysis: run the validated pipeline on every session in the panel.

Two independent questions, two nulls:

  1. SINGLE-CELL RHYTHMICITY  (ISI-shuffle null)
     Does a cell's autocorrelogram spectrum have excess power below 0.1 Hz?
     NOTE: shown by the bake-off to be non-specific -- it fires in the open
     field too, where grid/spatial tuning trivially produces slow rate
     modulation. Reported, but it is NOT sufficient evidence for the phenomenon.

  2. POPULATION SEQUENCES  (circular-shift null)  <-- the real test
     Do the cells ride a COMMON oscillation at staggered phases? The circular
     shift preserves every cell's own spectrum and destroys only the between-cell
     timing, so it tests exactly the paper's claim and nothing else.
"""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(__file__))
import loaders
import mec_oscillation as mo

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N_ISI = 100
N_SHIFT = 200
DECIM = 8          # 0.12 s -> ~0.96 s bins for the population test (signal is <0.1 Hz)

PANEL = [
    ("ctrl_wheel",  "POSITIVE CONTROL — original paper, wheel in darkness",
     lambda: loaders.load_ebrains("104638", "Wheel-HeadFixed")),
    ("ctrl_openfield", "CONTROL — original paper, open field (same units)",
     lambda: loaders.load_ebrains("104638", "OpenField")),
    ("d000053_mec", "DANDI 000053 — mouse MEC, VR straight track",
     lambda: loaders.load("000053", "sub-npI1/sub-npI1_ses-20190413_behavior+ecephys.nwb")),
    ("d001701_mec", "DANDI 001701 — mouse MEC, X-maze",
     lambda: loaders.load("001701", "sub-AppleBottom/sub-AppleBottom_ses-AppleBottom-DY05-g1_behavior+ecephys.nwb",
                          region_select="MEC")),
    ("d001701_v1", "DANDI 001701 — visual cortex, SAME session (region control)",
     lambda: loaders.load("001701", "sub-AppleBottom/sub-AppleBottom_ses-AppleBottom-DY05-g1_behavior+ecephys.nwb",
                          region_select="V1")),
    ("d000897_ec", "DANDI 000897 — MACAQUE entorhinal cortex, mental navigation",
     lambda: loaders.load("000897", "sub-amadeus/sub-amadeus_ses-08152019_behavior+ecephys.nwb")),
]


def analyze(tag, title, loader):
    print(f"\n{'=' * 88}\n{title}")
    sess = loader()
    print(f"  {sess.species} | {sess.region}")
    print(f"  {sess.n_units} units (>=0.5 Hz) | {sess.duration / 60:.1f} min")
    if sess.n_units < 10:
        print("  SKIP: <10 units")
        return None

    t, z, b = mo.build_activity(sess.spike_times, sess.t0, sess.t1)

    # 1. single-cell rhythmicity
    sc = mo.single_cell_test(sess.spike_times, sess.t0, sess.t1,
                             n_surrogates=N_ISI, progress=tqdm)
    frac = sc["sig"].mean()
    print(f"  [1] single-cell ultraslow rhythmicity (ISI-shuffle, FDR q<0.05): "
          f"{sc['sig'].sum()}/{sess.n_units} = {100 * frac:.0f}%")
    if sc["sig"].any():
        fp = sc["f_peak"][sc["sig"]]
        print(f"      median peak {np.median(fp):.4f} Hz  (period {1 / np.median(fp):.0f} s)")

    # 2. population sequences
    zz = z[::DECIM]
    pop = mo.population_sequence_test(zz, n_surrogates=N_SHIFT, progress=tqdm)
    print(f"  [2] POPULATION SEQUENCES (circular-shift null):")
    print(f"      PC1+PC2 variance : {100 * pop['pc12_var']:.1f}%  vs shuffled "
          f"{100 * pop['null_pc12_mean']:.1f}%   z={pop['pc12_z']:+.1f}  p={pop['p_pc12']:.4f}")
    print(f"      rotation index   : {pop['rotation']:.3f}  vs shuffled "
          f"{pop['null_rot_mean']:.3f}   z={pop['rot_z']:+.1f}  p={pop['p_rotation']:.4f}")

    order = np.argsort(-pop["theta"])
    ve = pop["var_exp"]

    # ---------------- figure ----------------
    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(3, 3, hspace=0.62, wspace=0.34,
                          left=0.07, right=0.96, top=0.87, bottom=0.07)
    fig.suptitle(f"{title}\n{sess.species} | {sess.region} | {sess.n_units} units | "
                 f"{sess.duration / 60:.0f} min", fontsize=12.5, y=0.965)

    ex = np.argsort(-sc["stat"])[:5]

    ax = fig.add_subplot(gs[0, 0])
    for k, i in enumerate(ex):
        ax.plot(t / 60, z[:, i] + k * 5, lw=0.7)
    ax.set_xlabel("time (min)"); ax.set_ylabel("z-scored rate (offset)")
    ax.set_title("(a) most-rhythmic cells\nsmoothed rate ($\\sigma$=5 s)", fontsize=10, pad=8)

    ax = fig.add_subplot(gs[0, 1])
    ml = int(300 / mo.BIN_SIZE)
    lags = np.arange(ml + 1) * mo.BIN_SIZE
    for i in ex:
        ax.plot(lags, mo.autocorr(b[:, i], ml), lw=0.9)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("lag (s)"); ax.set_ylabel("autocorrelation")
    ax.set_title("(b) autocorrelation", fontsize=10, pad=8)

    ax = fig.add_subplot(gs[0, 2])
    cnt = mo._bin_counts(sess.spike_times, sess.t0, sess.t1)
    for i in ex:
        f, p = mo.acg_spectrum(cnt[i])
        ax.semilogy(f, p, lw=0.9)
    ax.axvspan(*mo.BAND, color="grey", alpha=0.18)
    ax.set_xlim(0, 0.15)
    ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("PSD of autocorrelogram")
    ax.set_title("(c) autocorrelogram spectrum", fontsize=10, pad=8)

    ax = fig.add_subplot(gs[1, 0])
    bins = np.linspace(0, 0.1, 26)
    if (~sc["sig"]).any():
        ax.hist(sc["f_peak"][~sc["sig"]], bins=bins, color="lightgrey",
                label=f"n.s. (n={(~sc['sig']).sum()})")
    if sc["sig"].any():
        ax.hist(sc["f_peak"][sc["sig"]], bins=bins, color="tab:blue",
                label=f"rhythmic (n={sc['sig'].sum()})")
    ax.set_xlabel("peak frequency (Hz)"); ax.set_ylabel("# cells")
    ax.legend(fontsize=8)
    ax.set_title(f"(d) single-cell rhythmicity\n{100 * frac:.0f}% of cells", fontsize=10, pad=8)

    ax = fig.add_subplot(gs[1, 1])
    ax.plot(pop["scores"][:, 0], pop["scores"][:, 1], lw=0.4, color="k", alpha=0.6)
    ax.set_xlabel(f"PC1 ({100 * ve[0]:.0f}%)"); ax.set_ylabel(f"PC2 ({100 * ve[1]:.0f}%)")
    ax.set_title("(e) population trajectory", fontsize=10, pad=8)
    ax.set_aspect("equal", adjustable="datalim")

    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    verdict = ("SEQUENCES PRESENT" if (pop["p_pc12"] < 0.05 and pop["rot_z"] > 3)
               else "no population sequence")
    txt = (f"POPULATION TEST\n(circular-shift null, n={N_SHIFT})\n\n"
           f"PC1+PC2 variance\n"
           f"   data      {100 * pop['pc12_var']:.1f}%\n"
           f"   shuffled  {100 * pop['null_pc12_mean']:.1f}%\n"
           f"   z = {pop['pc12_z']:+.1f}   p = {pop['p_pc12']:.4f}\n\n"
           f"rotation index\n"
           f"   data      {pop['rotation']:.3f}\n"
           f"   shuffled  {pop['null_rot_mean']:.3f}\n"
           f"   z = {pop['rot_z']:+.1f}   p = {pop['p_rotation']:.4f}\n\n"
           f"→ {verdict}")
    ax.text(0.0, 0.98, txt, va="top", ha="left", fontsize=9.5, family="monospace",
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.6", fc="#eef3f8", ec="#7c9fc0"))

    ax = fig.add_subplot(gs[2, :])
    show = min(len(t), int(1200 / mo.BIN_SIZE))
    im = ax.imshow(z[:show, order].T, aspect="auto", cmap="viridis", vmin=-1, vmax=3,
                   origin="lower", extent=[t[0] / 60, t[show - 1] / 60, 0, sess.n_units])
    ax.set_xlabel("time (min)"); ax.set_ylabel("cell (sorted by PCA angle θ)")
    ax.set_title("(f) activity sorted by PC1/PC2 angle — repeating diagonal bands = periodic sequences",
                 fontsize=11, pad=8)
    plt.colorbar(im, ax=ax, label="z-scored rate", pad=0.01)

    out = f"{ROOT}/figures/{tag}.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"  -> {out}")

    summary = dict(
        tag=tag, title=title, species=sess.species, region=sess.region,
        n_units=int(sess.n_units), duration_min=float(sess.duration / 60),
        frac_rhythmic=float(frac), n_rhythmic=int(sc["sig"].sum()),
        median_peak_hz=float(np.median(sc["f_peak"][sc["sig"]])) if sc["sig"].any() else None,
        pc12_var=pop["pc12_var"], pc12_null=pop["null_pc12_mean"],
        pc12_z=pop["pc12_z"], p_pc12=pop["p_pc12"],
        rotation=pop["rotation"], rotation_null=pop["null_rot_mean"],
        rot_z=pop["rot_z"], p_rotation=pop["p_rotation"],
        sequences=bool(pop["p_pc12"] < 0.05 and pop["rot_z"] > 3),
    )
    with open(f"{ROOT}/results/{tag}.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary


if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    out = []
    for tag, title, loader in PANEL:
        if only and only != tag:
            continue
        r = analyze(tag, title, loader)
        if r:
            out.append(r)
    if out:
        print(f"\n\n{'=' * 88}\nSUMMARY")
        print(f"{'session':<46}{'units':>6}{'rhythmic':>10}{'PC12 z':>9}{'rot z':>8}  seq")
        print("-" * 88)
        for r in out:
            print(f"{r['tag']:<46}{r['n_units']:>6}{100 * r['frac_rhythmic']:>9.0f}%"
                  f"{r['pc12_z']:>9.1f}{r['rot_z']:>8.1f}  "
                  f"{'YES' if r['sequences'] else 'no'}")
