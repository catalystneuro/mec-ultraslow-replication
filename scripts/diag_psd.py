"""Diagnostic: compare the real PSD to the ISI-shuffle surrogate PSD, per cell.

The question: does the ultraslow band contain a *peak* that the surrogates lack,
or does smoothing simply put all power there for data and surrogates alike?
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

dsid, path = sys.argv[1], sys.argv[2]
n_surr = 30

sess = loaders.load(dsid, path)
print(f"{sess.n_units} units, {sess.duration/60:.1f} min")
t, z, b = mo.build_activity(sess.spike_times, sess.t0, sess.t1)
fs = 1 / mo.BIN_SIZE
rng = np.random.default_rng(0)

# surrogate PSDs
surr_psds = []
for k in range(n_surr):
    sur = [mo.isi_shuffle(st, sess.t0, sess.t1, rng) for st in sess.spike_times]
    _, zs, bs = mo.build_activity(sur, sess.t0, sess.t1)
    ps = []
    for i in range(sess.n_units):
        f, p = mo.compute_psd(bs[:, i], fs)
        ps.append(p)
    surr_psds.append(np.array(ps))
surr_psds = np.array(surr_psds)          # (n_surr, n_cells, n_freq)

data_psds = []
for i in range(sess.n_units):
    f, p = mo.compute_psd(b[:, i], fs)
    data_psds.append(p)
data_psds = np.array(data_psds)          # (n_cells, n_freq)

print("freq resolution:", f[1] - f[0], "Hz; max freq:", f[-1])

# ratio of data to mean surrogate, in band
msur = surr_psds.mean(axis=0)
ratio = data_psds / np.maximum(msur, 1e-30)
in_band = (f >= mo.BAND[0]) & (f <= mo.BAND[1])
peak_ratio = ratio[:, in_band].max(axis=1)
print(f"data/surrogate peak ratio in band: median {np.median(peak_ratio):.2f}, "
      f"max {peak_ratio.max():.2f}, n>2 = {(peak_ratio > 2).sum()}/{sess.n_units}")

# how much of the total power is in-band, for data vs surrogate?
tot_d = np.trapezoid(data_psds, f, axis=1)
bd = np.trapezoid(data_psds[:, in_band], f[in_band], axis=1) / tot_d
tot_s = np.trapezoid(msur, f, axis=1)
bs_ = np.trapezoid(msur[:, in_band], f[in_band], axis=1) / tot_s
print(f"fraction of power in ultraslow band: data {np.median(bd):.3f} vs surrogate {np.median(bs_):.3f}")
print("  ^ if these are both ~1.0, the band-power statistic cannot discriminate")

# plot the 6 cells with the biggest data/surrogate ratio
best = np.argsort(-peak_ratio)[:6]
fig, axes = plt.subplots(2, 3, figsize=(15, 7))
for ax, i in zip(axes.ravel(), best):
    ax.fill_between(f, np.percentile(surr_psds[:, i, :], 2.5, axis=0),
                    np.percentile(surr_psds[:, i, :], 97.5, axis=0),
                    color="grey", alpha=0.4, label="ISI-shuffle 95% CI")
    ax.plot(f, msur[i], color="k", lw=1, ls="--", label="surrogate mean")
    ax.plot(f, data_psds[i], color="tab:red", lw=1.2, label="data")
    ax.axvspan(*mo.BAND, color="tab:blue", alpha=0.08)
    ax.set_xlim(0, 0.12)
    ax.set_yscale("log")
    ax.set_title(f"unit {sess.unit_ids[i]}  peak ratio={peak_ratio[i]:.1f}", fontsize=9)
    ax.set_xlabel("Hz")
axes[0, 0].legend(fontsize=7)
fig.suptitle(f"{dsid} {os.path.basename(path)} — data vs ISI-shuffle PSD")
fig.tight_layout()
out = f"{ROOT}/figures/diag_psd_{dsid}.png"
fig.savefig(out, dpi=110)
print("->", out)
