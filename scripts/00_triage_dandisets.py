"""Triage candidate DANDI datasets with entorhinal cortex recordings.

Goal: find dandisets that could support a replication of Gonzalo Cogno et al.
(2024) minute-scale oscillatory sequences in MEC. Requirements:
  - spike times (units table) from MEC
  - long continuous recordings (>> 1 oscillation period; need >= ~10 min)
  - many simultaneously recorded units (need a population for sequences)
"""

from dandi.dandiapi import DandiAPIClient

CANDIDATES = [
    "000053",  # Mallory et al. 2021 - MEC linear track + open field (positive control from helpdesk)
    "000582",  # Conjunctive representation of position/direction/velocity in EC
    "000638",  # Hippocampus + EC dual region silicon probe
    "000897",  # Neupane/Fiete/Jazayeri - NHP entorhinal cortex, mental navigation
    "000943",  # Clark & Nolan 2024 - grid cells, path integration
    "001701",  # Aery Jones et al. 2026 - entorhinal cortex remote locations
]

client = DandiAPIClient()

for dsid in CANDIDATES:
    print("=" * 90)
    try:
        ds = client.get_dandiset(dsid, "draft")
        meta = ds.get_raw_metadata()
    except Exception as e:  # dandiset may be embargoed / unpublished
        print(f"{dsid}: could not fetch ({e})")
        continue

    print(f"{dsid}: {meta.get('name', '?')}")

    species = {s.get("name") for a in meta.get("about", []) for s in [a] if a.get("name")}
    # species is more reliably in assetsSummary
    summ = meta.get("assetsSummary", {})
    print(f"  species        : {[s.get('name') for s in summ.get('species', [])]}")
    print(f"  approach       : {[s.get('name') for s in summ.get('approach', [])]}")
    print(f"  measurement    : {[s.get('name') for s in summ.get('measurementTechnique', [])]}")
    print(f"  variables      : {summ.get('variableMeasured')}")
    print(f"  n subjects     : {summ.get('numberOfSubjects')}")
    print(f"  n files        : {summ.get('numberOfFiles')}")
    print(f"  total size     : {summ.get('numberOfBytes', 0) / 1e9:.1f} GB")

    # sample a few assets to see file naming / sizes
    n_show = 3
    for i, asset in enumerate(ds.get_assets()):
        if i >= n_show:
            break
        print(f"  asset[{i}]      : {asset.path}  ({asset.size / 1e6:.1f} MB)")
