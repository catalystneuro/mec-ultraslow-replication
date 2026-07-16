"""Is the 'oscillation' just the behavioural cycle?

The population sequence test flagged VISUAL CORTEX in the 001701 X-maze session
(rotation z=+5.8) more strongly than MEC (z=+0.7). V1 is supposed to be the
negative region control -- the original paper found no sequences there.

The suspicion: the animal runs laps, so the visual scene repeats periodically,
V1 cells tile the lap cycle, and the population traces a periodic sequence that
is behavioural, not intrinsic. The original study eliminated this by design
(wheel, darkness, no rewards, no landmarks).

Test: compute the ultraslow spectrum of the BEHAVIOUR itself with exactly the
same detector, and ask whether the neural population's peak frequency coincides
with the behavioural one.
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
from pynwb import NWBHDF5IO
from dandi.dandiapi import DandiAPIClient

import loaders
import mec_oscillation as mo

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASSET_1701 = "sub-AppleBottom/sub-AppleBottom_ses-AppleBottom-DY05-g1_behavior+ecephys.nwb"


def spectrum_of(x, fs):
    """Same detector as the neural data: spectrum of the z-scored autocorrelogram."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return mo.acg_spectrum(x - x.mean(), fs=fs)


def band_peak(f, p, band=mo.BAND):
    m = (f >= band[0]) & (f <= band[1])
    return f[m][np.argmax(p[m])]


def pop_pc1(sess):
    _, z, _ = mo.build_activity(sess.spike_times, sess.t0, sess.t1)
    x = z - z.mean(0, keepdims=True)
    c = np.cov(x, rowvar=False)
    w, v = np.linalg.eigh(c)
    pc1 = x @ v[:, np.argmax(w)]
    return pc1


# ---- behaviour from 001701 -------------------------------------------------
ds = DandiAPIClient().get_dandiset("001701", "draft")
asset = ds.get_asset_by_path(ASSET_1701)
f = remfile.File(asset.get_content_url(follow_redirects=1, strip_query=True),
                 disk_cache=remfile.DiskCache(loaders.CACHE))
io = NWBHDF5IO(file=h5py.File(f, "r"), load_namespaces=True)
nwb = io.read()

pos_iface = nwb.processing["behavior"]["Position"]
ss = list(pos_iface.spatial_series.values())[0]
pos = np.asarray(ss.data[:])
pos_t = (np.asarray(ss.timestamps[:]) if ss.timestamps is not None
         else ss.starting_time + np.arange(len(pos)) / ss.rate)
print(f"Position: {pos.shape}, t {pos_t[0]:.1f}-{pos_t[-1]:.1f} s")
io.close()

if pos.ndim > 1:
    pos = pos[:, 0]

# resample behaviour onto the 120 ms analysis grid
grid = np.arange(pos_t[0], pos_t[-1], mo.BIN_SIZE)
pos_r = np.interp(grid, pos_t, pos)
speed_r = np.abs(np.gradient(pos_r, mo.BIN_SIZE))
fs = 1 / mo.BIN_SIZE

fb, pb = spectrum_of(pos_r, fs)
fsp, psp = spectrum_of(speed_r, fs)
f_pos = band_peak(fb, pb)
f_spd = band_peak(fsp, psp)

print(f"\nBEHAVIOUR ultraslow peak:")
print(f"  position : {f_pos:.4f} Hz  (period {1/f_pos:.0f} s)")
print(f"  speed    : {f_spd:.4f} Hz  (period {1/f_spd:.0f} s)")

# ---- neural populations ----------------------------------------------------
print(f"\nNEURAL population (PC1) ultraslow peak:")
rows = []
for reg in ("MEC", "V1"):
    sess = loaders.load("001701", ASSET_1701, region_select=reg)
    pc1 = pop_pc1(sess)
    fn, pn = spectrum_of(pc1, fs)
    f_pc1 = band_peak(fn, pn)
    print(f"  {reg:<4} (n={sess.n_units:>3}) : {f_pc1:.4f} Hz  (period {1/f_pc1:.0f} s)")
    rows.append((reg, sess.n_units, fn, pn, f_pc1))

# ---- figure ----------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
ax.plot(grid[:5000] / 60, pos_r[:5000], color="k", lw=0.8)
ax.set_xlabel("time (min)"); ax.set_ylabel("position")
ax.set_title("001701 X-maze: the animal runs a repeating lap cycle", fontsize=11)

ax = axes[1]
ax.semilogy(fb, pb / pb.max(), color="k", lw=1.6, label=f"behaviour: position ({f_pos:.3f} Hz)")
for (reg, n, fn, pn, f_pc1), col in zip(rows, ["tab:red", "tab:blue"]):
    ax.semilogy(fn, pn / pn.max(), color=col, lw=1.2,
                label=f"{reg} population PC1 ({f_pc1:.3f} Hz)")
ax.axvline(f_pos, color="k", ls=":", lw=1)
ax.axvspan(*mo.BAND, color="grey", alpha=0.12)
ax.set_xlim(0, 0.12)
ax.set_xlabel("frequency (Hz)"); ax.set_ylabel("normalized PSD")
ax.legend(fontsize=9)
ax.set_title("Population 'oscillation' sits at the behavioural lap frequency", fontsize=11)

fig.tight_layout()
out = f"{ROOT}/figures/behavior_control_001701.png"
fig.savefig(out, dpi=120)
print(f"\n-> {out}")
