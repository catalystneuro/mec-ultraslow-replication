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


# 000690 sessions with medial entorhinal coverage, from a scan of all 25 sessions'
# electrodes tables (scan_000690_regions.py). All 11 also carry visual cortex on the
# same probes; only 4 carry parasubiculum, and thinly.
OPENSCOPE_MEC_SESSIONS = [
    "sub-695762/sub-695762_ses-1317448357.nwb",   # ENTm=348 PAR=77  VIS=304
    "sub-695763/sub-695763_ses-1317661297.nwb",   # ENTm=268 PAR=0   VIS=363
    "sub-695764/sub-695764_ses-1311204385.nwb",   # ENTm=266 PAR=0   VIS=348
    "sub-699321/sub-699321_ses-1312636156.nwb",   # ENTm=379 PAR=10  VIS=256
    "sub-699322/sub-699322_ses-1317198704.nwb",   # ENTm=57  PAR=0   VIS=207
    "sub-702134/sub-702134_ses-1324803287.nwb",   # ENTm=411 PAR=15  VIS=329
    "sub-714612/sub-714612_ses-1325995398.nwb",   # ENTm=360 PAR=10  VIS=289
    "sub-714615/sub-714615_ses-1325748772.nwb",   # ENTm=114 PAR=0   VIS=260
    "sub-716465/sub-716465_ses-1330689294.nwb",   # ENTm=99  PAR=0   VIS=276
    "sub-717437/sub-717437_ses-1333971345.nwb",   # ENTm=250 PAR=0   VIS=352
    "sub-719667/sub-719667_ses-1333741475.nwb",   # ENTm=343 PAR=33  VIS=298
]


def load_000690(asset_path, region_select="MEC", min_rate=0.5, trim_start_s=300.0,
                quality=("good",)):
    """Allen OpenScope Vision2Hippocampus: the best MEC coverage on DANDI.

    Why it is in the panel. A single session carries up to 726 good MEC units --
    nearly double the paper's own Neuropixels sessions (404 and 273) -- AND, on the
    same probes at the same moment, the paper's two negative-control regions:
    parasubiculum (PAR) and visual cortex (VIS). It therefore reproduces the design
    of the paper's Fig. 5 (MEC 15/27, PaS 0/25, VIS 0/19) WITHIN a session, whereas
    the original compared across different mice and days.

    That within-session structure is what makes it worth having despite the
    condition being wrong. The recurring failure mode in this panel is that apparent
    sequences turn out to be behavioural or task structure (V1's lap cycle in 001701,
    task-engagement blocks in the macaque). Here the stimulus is identical across
    simultaneously recorded regions by construction, so a MEC-vs-VIS contrast cannot
    be manufactured by stimulus block structure: if it could, VIS would show it at
    least as strongly.

    What it CANNOT do: test the paper's condition. There is no darkness and no
    task-free block -- stimuli tile t=247-7488 s with no gap >30 s. So an absolute
    result here needs the same skepticism applied to V1 and the macaque, and a null
    is ambiguous between "sequences reset under sensory drive" (the authors' own
    speculation) and "no effect". The RELATIVE, within-session region contrast is
    the robust part.

    Units map to regions via peak_channel_id -> electrodes.id -> location, using
    Allen CCF acronyms (ENTm* = medial entorhinal, PAR = parasubiculum, VIS* = any
    visual area). Only 'good' units are kept, matching the EBRAINS control loader.
    """
    io, nwb = _open("000690", asset_path)

    el = nwb.electrodes.to_dataframe()
    loc_by_id = {int(i): str(l) for i, l in zip(el.index, el["location"])}

    def is_region(s):
        s = s.strip()
        if region_select == "MEC":
            return s.startswith("ENTm")          # medial only; excludes ENTl
        if region_select == "PaS":
            return s.startswith("PAR")
        if region_select == "V1":
            return s.startswith("VIS")
        return False

    units = nwb.units
    pk = np.asarray(units["peak_channel_id"][:])
    ulocs = np.array([loc_by_id.get(int(x), "?") for x in pk])
    qual = (np.asarray(units["quality"][:], dtype=str)
            if "quality" in units.colnames else np.array(["good"] * len(units)))

    keep = [i for i in range(len(units))
            if is_region(ulocs[i]) and qual[i] in quality]

    spike_times = [np.asarray(units["spike_times"][i]) for i in keep]
    unit_ids = np.array(keep)
    if not spike_times:
        io.close()
        raise ValueError(f"no {region_select} units in {asset_path}")

    t0, t1 = _span(spike_times)
    t0 += trim_start_s
    spike_times = [s[(s >= t0) & (s <= t1)] for s in spike_times]
    spike_times, unit_ids = _filter_rate(spike_times, unit_ids, t0, t1, min_rate)

    subject = nwb.subject.subject_id if nwb.subject else "?"
    io.close()
    region = {"MEC": "MEC (ENTm, Allen CCF)",
              "PaS": "Parasubiculum (PAR)",
              "V1": "Visual cortex (VIS*)"}[region_select]
    return Session(dandiset="000690", asset_path=asset_path, subject=subject,
                   species="Mus musculus", region=region,
                   spike_times=spike_times, unit_ids=unit_ids, t0=t0, t1=t1,
                   description="passive visual stimulation, head-fixed (no darkness)")


def load_000690_arousal(asset_path):
    """Pupil area and running speed for a 000690 session, from its `_image.nwb`.

    This is the control the original paper concedes it cannot run. The paper
    attributes the single-cell oscillation to "ascending neuromodulatory
    arousal-associated brain-stem circuits" and cites infraslow PUPIL oscillation
    (Blasiak 2013) and vascular contributions (Drew 2020) -- but has no pupil,
    vascular or arousal measure anywhere, and 0.006-0.008 Hz sits squarely in the
    arousal/vasomotion range.

    Returns (t_pupil, pupil_area, t_run, running_speed); any element may be None.
    """
    img = asset_path.replace(".nwb", "_image.nwb")
    io, nwb = _open("000690", img)
    tp = pupil = tr = run = None

    eye = nwb.acquisition.get("EyeTracking")
    if eye is not None:
        # EllipseSeries: prefer the pupil ellipse's area
        ts = None
        for name in ("pupil_tracking", "eye_tracking"):
            ts = getattr(eye, name, None) or (eye.spatial_series.get(name)
                                              if hasattr(eye, "spatial_series") else None)
            if ts is not None:
                break
        if ts is not None:
            area = getattr(ts, "area", None)
            pupil = np.asarray(area[:]) if area is not None else np.asarray(ts.data[:])
            if pupil.ndim > 1:
                pupil = pupil[:, 0]
            tp = np.asarray(ts.timestamps[:]) if ts.timestamps is not None else None

    rp = nwb.processing.get("running")
    if rp is not None and "running_speed" in rp.data_interfaces:
        rs = rp["running_speed"]
        run = np.asarray(rs.data[:]).squeeze()
        tr = (np.asarray(rs.timestamps[:]) if rs.timestamps is not None
              else rs.starting_time + np.arange(len(run)) / rs.rate)

    io.close()
    return tp, pupil, tr, run


def load_000552(asset_path, state="Awake", min_rate=0.5, trim_start_s=300.0):
    """Hippocampal CA1 during long TASK-FREE rest (Huszar et al., dandiset 000552).

    Why this dataset is here: every entorhinal dataset on DANDI runs a task from
    start to finish, so none can be asked the paper's question in the paper's
    condition (no task, no scheduled rewards, minimal sensory drive). 000552 is
    the archive's best available approximation -- multi-hour sessions in which
    PRE/POST home-cage blocks bracket a maze epoch, with 100-376 units. It trades
    the region (CA1, not MEC) for the condition.

    That trade makes this a REGION control in the paper's own logic, not a
    replication: the paper's contrast is that MEC carries a low-dimensional code
    where ~94% of cells lock to the sequence, whereas hippocampus carries a
    high-dimensional code in which reported sequences involve ~5% of the network.
    A CA1 null is therefore the EXPECTED outcome and corroborates rather than
    tests. It does, however, bear on the authors' open question about sleep and
    rest, and it is the only public data that can be asked in a task-free state.

    Epoch selection: the epochs table is UNLABELLED (only start/stop), so
    task-free epochs are identified empirically as those containing no trials.
    Brain state comes from processing/behavior/SleepStates. We take the single
    longest continuous bout of `state` inside a task-free epoch -- concatenating
    bouts would destroy the phase continuity a minute-scale sequence lives in.
    As in the ephys pipeline elsewhere, the first 5 min are dropped for arousal.

    Caveat: electrodes.location is 'unknown' throughout; CA1 comes from the
    publication, not the file. Darkness is undocumented (home-cage rest is
    presumably dim/dark, but no metadata asserts it).
    """
    io, nwb = _open("000552", asset_path)

    ep = nwb.epochs.to_dataframe()
    trials = nwb.trials.to_dataframe() if nwb.trials is not None else None
    ts = np.asarray(trials["start_time"]) if trials is not None else np.array([])
    ts = ts[np.isfinite(ts)]

    task_free = [(float(r.start_time), float(r.stop_time)) for r in ep.itertuples()
                 if not ((ts >= r.start_time) & (ts < r.stop_time)).any()]

    states = nwb.processing["behavior"]["SleepStates"].to_dataframe()
    best = (0.0, None, None)
    for r in states.itertuples():
        if str(r.label) != state:
            continue
        for a0, b0 in task_free:
            x, y = max(float(r.start_time), a0), min(float(r.stop_time), b0)
            if y - x > best[0]:
                best = (y - x, x, y)
    if best[1] is None:
        io.close()
        raise ValueError(f"no task-free '{state}' bout in {asset_path}")

    t0, t1 = best[1] + trim_start_s, best[2]

    units = nwb.units
    spike_times, unit_ids = [], []
    for i in range(len(units)):
        s = np.asarray(units["spike_times"][i])
        s = s[(s >= t0) & (s <= t1)]
        spike_times.append(s)
        unit_ids.append(i)
    unit_ids = np.array(unit_ids)
    spike_times, unit_ids = _filter_rate(spike_times, unit_ids, t0, t1, min_rate)

    subject = nwb.subject.subject_id if nwb.subject else "?"
    io.close()
    return Session(dandiset="000552", asset_path=asset_path, subject=subject,
                   species="Mus musculus",
                   region="CA1 (from publication; electrodes.location unannotated)",
                   spike_times=spike_times, unit_ids=unit_ids, t0=t0, t1=t1,
                   description=f"task-free {state} rest ({(t1-t0)/60:.0f} min)")


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
