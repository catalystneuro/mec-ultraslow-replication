"""Stream the units table + session metadata from candidate sessions.

We only touch small metadata datasets, so this is cheap even for 30-70 GB assets.
Reports what the replication needs: recording duration, number of simultaneously
recorded units, brain regions, and available behavioral streams.
"""

import os
import sys
import warnings

import h5py
import numpy as np
import remfile
from dandi.dandiapi import DandiAPIClient
from pynwb import NWBHDF5IO

warnings.filterwarnings("ignore")

CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cache")


def open_asset(asset):
    disk_cache = remfile.DiskCache(CACHE)
    f = remfile.File(asset.get_content_url(follow_redirects=1, strip_query=True), disk_cache=disk_cache)
    io = NWBHDF5IO(file=h5py.File(f, "r"), load_namespaces=True)
    return io, io.read()


def inspect(dsid, n_sessions=3, path_filter=None):
    client = DandiAPIClient()
    ds = client.get_dandiset(dsid, "draft")

    assets = [a for a in ds.get_assets() if a.path.endswith(".nwb")]
    if path_filter:
        assets = [a for a in assets if path_filter in a.path]
    assets = assets[:n_sessions]

    print("#" * 95)
    print(f"# DANDISET {dsid}")
    print("#" * 95)

    for asset in assets:
        print(f"\n--- {asset.path}  ({asset.size / 1e6:.0f} MB)")
        io, nwb = open_asset(asset)

        print(f"    session_description: {str(nwb.session_description)[:100]}")
        print(f"    subject            : {nwb.subject.subject_id if nwb.subject else '?'} "
              f"species={nwb.subject.species if nwb.subject else '?'}")

        if nwb.units is None:
            print("    UNITS: none")
            io.close()
            continue

        n_units = len(nwb.units)
        st = nwb.units["spike_times"]
        # session span from spike times (cheap: read first/last per unit is not possible,
        # so read all spike times -- units tables are small relative to raw data)
        all_first, all_last, counts = [], [], []
        for i in range(n_units):
            s = st[i]
            counts.append(len(s))
            if len(s):
                all_first.append(s[0])
                all_last.append(s[-1])
        t0, t1 = float(np.min(all_first)), float(np.max(all_last))
        dur = t1 - t0
        rates = np.array(counts) / max(dur, 1e-9)

        print(f"    UNITS              : {n_units}")
        print(f"    spike-time span    : {t0:.1f} -> {t1:.1f} s  ({dur / 60:.1f} min)")
        print(f"    firing rate        : median {np.median(rates):.2f} Hz, "
              f"n>0.1Hz = {(rates > 0.1).sum()}, n>0.5Hz = {(rates > 0.5).sum()}")
        print(f"    units columns      : {list(nwb.units.colnames)}")

        # brain region info
        for col in ("location", "brain_region", "region", "area"):
            if col in nwb.units.colnames:
                vals, cnt = np.unique(np.asarray(nwb.units[col][:], dtype=str), return_counts=True)
                print(f"    units.{col}       : {dict(zip(vals, cnt))}")
        if nwb.electrodes is not None and "location" in nwb.electrodes.colnames:
            vals, cnt = np.unique(np.asarray(nwb.electrodes["location"][:], dtype=str), return_counts=True)
            print(f"    electrode locations: {dict(zip(vals, cnt))}")

        print(f"    acquisition        : {list(nwb.acquisition.keys())}")
        print(f"    processing         : { {k: list(v.data_interfaces.keys()) for k, v in nwb.processing.items()} }")
        print(f"    intervals          : {list(nwb.intervals.keys()) if nwb.intervals else []}")
        if nwb.trials is not None:
            print(f"    trials             : {len(nwb.trials)} | cols={list(nwb.trials.colnames)}")
        io.close()


if __name__ == "__main__":
    dsid = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    pf = sys.argv[3] if len(sys.argv) > 3 else None
    inspect(dsid, n, pf)
