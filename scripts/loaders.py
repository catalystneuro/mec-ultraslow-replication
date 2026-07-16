"""Session loaders for the three dandisets in the replication panel.

Each returns a common Session record so the analysis pipeline is dataset-agnostic.
"""

import os
import warnings
from dataclasses import dataclass, field

import h5py
import numpy as np
import remfile
from dandi.dandiapi import DandiAPIClient
from pynwb import NWBHDF5IO

warnings.filterwarnings("ignore")

CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cache")

# Which electrode 'location' strings count as medial entorhinal cortex
MEC_PATTERNS = ("entorhinal area medial", "mec", "entorhinal")


@dataclass
class Session:
    dandiset: str
    asset_path: str
    subject: str
    species: str
    region: str
    spike_times: list = field(repr=False)
    unit_ids: np.ndarray = field(repr=False)
    t0: float = 0.0
    t1: float = 0.0
    speed: tuple = None       # (t, v) if available
    description: str = ""

    @property
    def n_units(self):
        return len(self.spike_times)

    @property
    def duration(self):
        return self.t1 - self.t0


EBRAINS_URLS = {
    "104638": "https://rgw.cscs.ch/ebrains:d-a25f722e-f9b1-4f32-80bc-37206c2d1536/data/ecephys/104638/104638_2023-Feb-10_11-04-33.nwb",
    "102335": "https://rgw.cscs.ch/ebrains:d-a25f722e-f9b1-4f32-80bc-37206c2d1536/data/ecephys/102335/102335_2022-Jul-01_09-28-53.nwb",
}


def load_ebrains(mouse, task, min_rate=0.5, trim_start_s=300.0, quality=("good",)):
    """The ORIGINAL Gonzalo Cogno et al. Neuropixels data from EBRAINS
    (doi:10.25493/SKKX-4W3).

    Each session contains two epochs recorded from the same MEC units:
      task='Wheel-HeadFixed' -> rotating wheel in darkness -> effect KNOWN PRESENT.
                                This is the positive control: it is how we prove the
                                detector can find the effect before trusting any null.
      task='OpenField'       -> free foraging -> effect UNKNOWN. The paper's Methods
                                state these trials "were not used in the present
                                study", so this is an untested condition, NOT a
                                validated negative control. Results on it are new
                                observations bearing on the authors' open question
                                about free exploration.

    Following the reference implementation, we drop the first 5 min of the epoch and
    (for mouse 104638) keep only 'good' units.
    """
    f = remfile.File(EBRAINS_URLS[mouse], disk_cache=remfile.DiskCache(CACHE))
    io = NWBHDF5IO(file=h5py.File(f, "r"), load_namespaces=True)
    nwb = io.read()

    ep = nwb.epochs.to_dataframe()
    row = ep[ep["task"] == task].iloc[0]
    t0, t1 = float(row["start_time"]) + trim_start_s, float(row["stop_time"])

    units = nwb.units
    qual = np.asarray(units["cluster_quality"][:], dtype=str)
    # Mouse 102335 carries TWO probes (MEC + hippocampal CA1); 104638 carries one
    # (MEC). Restrict to the MEC probe via each unit's electrode-group location,
    # or the CA1 units contaminate the MEC population and wash out the sequence.
    uloc = [units["electrode_group"][i].location for i in range(len(units))]
    keep = [i for i in range(len(units))
            if qual[i] in quality and "MEC" in uloc[i]]

    spike_times, unit_ids = [], []
    for i in keep:
        s = np.asarray(units["spike_times"][i])
        s = s[(s >= t0) & (s <= t1)]
        spike_times.append(s)
        unit_ids.append(i)
    unit_ids = np.array(unit_ids)

    spike_times, unit_ids = _filter_rate(spike_times, unit_ids, t0, t1, min_rate)
    io.close()

    return Session(dandiset="EBRAINS (original paper)", asset_path=f"{mouse} / {task}",
                   subject=mouse, species="Mus musculus",
                   region="MEC - Layer 2, 3, 5 (Neuropixels 2.0)",
                   spike_times=spike_times, unit_ids=unit_ids, t0=t0, t1=t1,
                   description=task)


def _open(dsid, asset_path):
    ds = DandiAPIClient().get_dandiset(dsid, "draft")
    asset = ds.get_asset_by_path(asset_path)
    f = remfile.File(asset.get_content_url(follow_redirects=1, strip_query=True),
                     disk_cache=remfile.DiskCache(CACHE))
    io = NWBHDF5IO(file=h5py.File(f, "r"), load_namespaces=True)
    return io, io.read()


def list_assets(dsid, contains=None, endswith=".nwb"):
    ds = DandiAPIClient().get_dandiset(dsid, "draft")
    paths = [a.path for a in ds.get_assets() if a.path.endswith(endswith)]
    if contains:
        paths = [p for p in paths if contains in p]
    return sorted(paths)


def _span(spike_times):
    firsts = [s[0] for s in spike_times if len(s)]
    lasts = [s[-1] for s in spike_times if len(s)]
    return float(min(firsts)), float(max(lasts))


def _filter_rate(spike_times, unit_ids, t0, t1, min_rate):
    dur = t1 - t0
    keep = [i for i, s in enumerate(spike_times) if len(s) / dur >= min_rate]
    return [spike_times[i] for i in keep], unit_ids[keep]


def load(dsid, asset_path, min_rate=0.5, region_select="MEC"):
    """Load one session, keeping only units from `region_select` above min_rate (Hz)."""
    io, nwb = _open(dsid, asset_path)

    subject = nwb.subject.subject_id if nwb.subject else "?"
    species = str(nwb.subject.species) if nwb.subject else "?"
    desc = str(nwb.session_description)[:80]
    speed = None

    if dsid == "000897":
        # units live in processing/ecephys/units (non-standard location).
        # The whole probe targets entorhinal cortex (dandiset title), and
        # electrodes.location is 'unknown', so we take all units.
        units = nwb.processing["ecephys"]["units"]
        spike_times = [np.asarray(units["spike_times"][i]) for i in range(len(units))]
        unit_ids = np.arange(len(units))
        region = "EC (all probe units; electrode location not annotated)"

    else:
        units = nwb.units
        spike_times = [np.asarray(units["spike_times"][i]) for i in range(len(units))]
        unit_ids = np.arange(len(units))

        if dsid == "001701":
            # map each unit -> its electrode -> brain region.
            # This dandiset records MEC, parasubiculum AND visual cortex on the
            # same probe, reproducing the paper's own regional control (effect
            # specific to MEC, absent in visual cortex). CRITICAL: two location
            # naming conventions coexist across subjects -- spelled-out ("Entorhinal
            # area medial part dorsal zone", "Primary visual area") for some, and
            # Allen atlas abbreviations ("ENTm2/ENTm3/ENTm5", "VISpl5") for others.
            # Matching only the spelled-out form silently drops ~half the sessions
            # (subjects Curie, Lamarr, ...), which all use the ENTm codes.
            def is_region(loc):
                s = loc.lower().strip()
                if region_select == "MEC":
                    # ENTm = entorhinal area, medial part (Allen). Exclude ENTl
                    # (lateral). "entm" and the spelled-out MEC both qualify.
                    return ("entorhinal area medial" in s
                            or s.startswith("entm") or s == "mec"
                            or s.startswith("mec"))
                if region_select == "PaS":
                    return "parasubiculum" in s or s.startswith("pas")
                if region_select == "V1":
                    # any visual cortical area (VISp primary, VISpl posterolateral, ...)
                    return "visual area" in s or s.startswith("vis")
                return False

            # Read each unit's location from its electrode sub-table DIRECTLY.
            # Do NOT index a positional location array by unit_elec[i].index --
            # those are electrode *ids*, which are non-contiguous on multi-shank
            # probes (e.g. 0-95, 1280-1375, ...), so positional indexing silently
            # returns the WRONG electrodes' locations.
            unit_elec = units["electrodes"]
            keep = []
            for i in range(len(units)):
                sub = unit_elec[i]
                locs = [str(x) for x in sub["location"].values] if len(sub) else []
                if any(is_region(l) for l in locs):
                    keep.append(i)
            spike_times = [spike_times[i] for i in keep]
            unit_ids = unit_ids[np.array(keep, dtype=int)] if keep else unit_ids[:0]
            region = {"MEC": "MEC (Entorhinal area medial part, dorsal zone)",
                      "PaS": "Parasubiculum",
                      "V1": "Visual cortex (V1 + posterolateral)"}[region_select]

        elif dsid == "000053":
            # Neuropixels probe targeted at MEC; electrodes.location is 'unknown'
            region = "MEC (Neuropixels, targeted; location not annotated)"

        else:
            region = "?"

        # behavioral speed, where available
        beh = nwb.processing.get("behavior")
        if beh is not None:
            for key in ("body_speed", "speed"):
                if key in beh.data_interfaces:
                    ts = beh[key]
                    ts = ts if not hasattr(ts, "time_series") else list(ts.time_series.values())[0]
                    v = np.asarray(ts.data[:]).squeeze()
                    if ts.timestamps is not None:
                        tt = np.asarray(ts.timestamps[:])
                    else:
                        tt = ts.starting_time + np.arange(len(v)) / ts.rate
                    speed = (tt, v)
                    break

    t0, t1 = _span(spike_times)
    spike_times, unit_ids = _filter_rate(spike_times, unit_ids, t0, t1, min_rate)

    io.close()
    return Session(dandiset=dsid, asset_path=asset_path, subject=subject,
                   species=species, region=region, spike_times=spike_times,
                   unit_ids=unit_ids, t0=t0, t1=t1, speed=speed, description=desc)
