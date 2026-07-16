"""Fast metadata scan: which 001701 sessions actually have MEC / V1 units?

Reads only the units' electrode-location column (cheap), so we know the true
denominator before committing to a full re-analysis. Writes a manifest.
"""
import json
import os
import sys
import warnings
from collections import Counter

import h5py
import numpy as np
import remfile
from pynwb import NWBHDF5IO

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import loaders

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def is_mec(s):
    s = s.lower().strip()
    return "entorhinal area medial" in s or s.startswith("entm") or s.startswith("mec")


def is_v1(s):
    s = s.lower().strip()
    return "visual area" in s or s.startswith("vis")


ds = loaders.DandiAPIClient().get_dandiset("001701", "draft")
paths = loaders.list_assets("001701")
print(f"scanning {len(paths)} sessions...")

manifest = {}
subj_mec = Counter()
for i, p in enumerate(paths):
    a = ds.get_asset_by_path(p)
    rf = remfile.File(a.get_content_url(follow_redirects=1, strip_query=True),
                      disk_cache=remfile.DiskCache(loaders.CACHE))
    io = NWBHDF5IO(file=h5py.File(rf, "r"), load_namespaces=True)
    nwb = io.read()
    u = nwb.units
    ue = u["electrodes"]
    locs = []
    for j in range(len(u)):
        sub = ue[j]
        locs.append(str(sub["location"].values[0]) if len(sub) else "none")
    n_mec = sum(is_mec(l) for l in locs)
    n_v1 = sum(is_v1(l) for l in locs)
    subj = p.split("/")[0]
    manifest[p] = dict(subject=subj, n_units=len(u), n_mec=n_mec, n_v1=n_v1)
    if n_mec >= 20:
        subj_mec[subj] += 1
    io.close()
    if (i + 1) % 25 == 0:
        print(f"  {i+1}/{len(paths)}", flush=True)

with open(f"{ROOT}/results/manifest_001701.json", "w") as fh:
    json.dump(manifest, fh, indent=2)

n_mec_sessions = sum(1 for v in manifest.values() if v["n_mec"] >= 20)
n_v1_sessions = sum(1 for v in manifest.values() if v["n_v1"] >= 20)
print(f"\nsessions with >=20 MEC units: {n_mec_sessions}/{len(paths)}")
print(f"sessions with >=20 V1 units : {n_v1_sessions}/{len(paths)}")
print(f"subjects contributing MEC   : {dict(subj_mec)}")
