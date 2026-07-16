"""Is the macaque EC 'oscillation' intrinsic, or task/arousal structure?

8/15 macaque entorhinal sessions pass the sequence test (some at p~1e-21). Before
calling that a cross-species replication, we must rule out the confound that
sank the V1 result: task structure imposing periodicity. The mouse paper's whole
claim is that the sequences are INTRINSIC -- they "transcended epochs of
immobility" and did not track behaviour.

Macaque sessions are 2.5-5.7 h of a trial-structured mental-navigation task. Two
diagnostics per session:

  1. SPECTRAL COINCIDENCE. Compute the population PC1 oscillation spectrum and the
     spectra of task regressors (trial-onset density, target value, success).
     Same detector on each. If the population peak sits on a task-regressor peak,
     the 'oscillation' is task structure.

  2. DOES IT SURVIVE THE GAPS. The task has long inter-trial gaps. If the
     sequences are intrinsic they should persist through gaps (as in the mouse,
     through immobility). If they only appear during dense task engagement, they
     are task-driven.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import h5py
import remfile
from dandi.dandiapi import DandiAPIClient
from pynwb import NWBHDF5IO

import loaders
import mec_oscillation as mo

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# strongest sessions from the run
SESSIONS = [
    ("sub-mahler/sub-mahler_ses-04122021_behavior+ecephys.nwb", "mahler 0412 (p~1e-21)"),
    ("sub-mahler/sub-mahler_ses-04142021_behavior+ecephys.nwb", "mahler 0414"),
    ("sub-amadeus/sub-amadeus_ses-08292019_behavior+ecephys.nwb", "amadeus 0829"),
]


def spec(x, fs):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 100 or x.std() == 0:
        return np.array([0.0]), np.array([0.0])
    return mo.acg_spectrum(x - x.mean(), fs=fs)


def band_peak(f, p):
    m = (f >= mo.BAND[0]) & (f <= mo.BAND[1])
    return f[m][np.argmax(p[m])] if m.any() and p[m].max() > 0 else np.nan


def pop_pc1(z):
    x = z - z.mean(0, keepdims=True)
    c = np.cov(x, rowvar=False)
    w, v = np.linalg.eigh(c)
    return x @ v[:, np.argmax(w)]


client = DandiAPIClient()
ds = client.get_dandiset("000897", "draft")

fig, axes = plt.subplots(len(SESSIONS), 2, figsize=(15, 4.2 * len(SESSIONS)))
if len(SESSIONS) == 1:
    axes = axes[None, :]

for row, (path, label) in enumerate(SESSIONS):
    sess = loaders.load("000897", path)
    t, z, b = mo.build_activity(sess.spike_times, sess.t0, sess.t1)
    fs = 1.0 / mo.BIN_SIZE
    pc1 = pop_pc1(z)
    fp, pp = spec(pc1, fs)
    f_pop = band_peak(fp, pp)

    # task regressors on the same time grid
    asset = ds.get_asset_by_path(path)
    f = remfile.File(asset.get_content_url(follow_redirects=1, strip_query=True),
                     disk_cache=remfile.DiskCache(loaders.CACHE))
    io = NWBHDF5IO(file=h5py.File(f, "r"), load_namespaces=True)
    nwb = io.read()
    tr = nwb.trials.to_dataframe()
    io.close()

    starts = tr["start_time"].values
    starts = starts[np.isfinite(starts)]
    grid = t
    # trial-onset density: count starts per bin, then smooth like the neural data
    onset = np.histogram(starts, bins=np.append(grid, grid[-1] + mo.BIN_SIZE))[0].astype(float)
    from scipy.ndimage import gaussian_filter1d
    onset_s = gaussian_filter1d(onset, mo.SIGMA / mo.BIN_SIZE)
    fo, po = spec(onset_s, fs)
    f_task = band_peak(fo, po)

    # target value over time (piecewise-constant across each trial)
    target = np.full(len(grid), np.nan)
    if "target" in tr.columns:
        for s, e, tg in zip(tr["start_time"], tr["stop_time"], tr["target"]):
            if np.isfinite(s) and np.isfinite(e):
                target[(grid >= s) & (grid < e)] = tg
    ftg, ptg = spec(np.nan_to_num(target, nan=np.nanmean(target)), fs)
    f_target = band_peak(ftg, ptg)

    print(f"{label:<26} pop peak={f_pop:.4f}Hz ({1/f_pop:.0f}s) | "
          f"trial-onset peak={f_task:.4f}Hz | target peak={f_target:.4f}Hz")

    # left: PC1 trace + trial onset density
    ax = axes[row, 0]
    tt = (grid - grid[0]) / 60
    ax.plot(tt, (pc1 - pc1.mean()) / pc1.std(), color="tab:blue", lw=0.5, label="population PC1")
    ax.plot(tt, (onset_s - onset_s.mean()) / (onset_s.std() + 1e-9) - 6,
            color="k", lw=0.6, label="trial-onset density")
    ax.set_xlim(0, min(60, tt[-1]))
    ax.set_xlabel("time (min)"); ax.set_yticks([])
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title(f"{label}: first hour", fontsize=10)

    # right: spectra overlaid
    ax = axes[row, 1]
    for (fx, px, c, lab, fpk) in [(fp, pp, "tab:blue", "population PC1", f_pop),
                                  (fo, po, "k", "trial-onset density", f_task),
                                  (ftg, ptg, "tab:orange", "target value", f_target)]:
        if px.max() > 0:
            ax.semilogy(fx, px / px.max(), color=c, lw=1.2, label=f"{lab} ({fpk:.3f} Hz)")
    ax.axvspan(*mo.BAND, color="grey", alpha=0.1)
    ax.set_xlim(0, 0.12)
    ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("norm. PSD")
    ax.legend(fontsize=8)
    ax.set_title(f"{label}: population vs task spectra", fontsize=10)

fig.suptitle("Macaque EC: is the population oscillation intrinsic or task-driven?", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = f"{ROOT}/figures/macaque_task_control.png"
fig.savefig(out, dpi=115)
print("->", out)
