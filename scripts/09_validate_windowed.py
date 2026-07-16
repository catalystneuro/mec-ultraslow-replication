"""Validate the windowed sequence test on the paper's OWN two Neuropixels mice.

The paper reports periodic sequences in BOTH Neuropixels wheel sessions
(mouse 104638, Fig. 2f; mouse 102335, Extended Data Fig. 4g). A usable detector
must fire on BOTH. The whole-session test only fired on 104638 -- it returned
rot_z = -5.2 on 102335, a false negative on known-positive data, because the
sequences are intermittent within the session.

Also note: the reference implementation uses good+MUA units for 102335 and
good-only for 104638, matching the paper's own unit counts (410 and 469).
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import loaders
import mec_oscillation as mo

DECIM = 8
BIN_S = mo.BIN_SIZE * DECIM

CASES = [
    ("104638", "Wheel-HeadFixed", ("good",), "POSITIVE (paper Fig. 2f)"),
    ("102335", "Wheel-HeadFixed", ("good", "mua"), "POSITIVE (paper Ext. Fig. 4g)"),
    ("104638", "OpenField", ("good",), "untested condition"),
    ("102335", "OpenField", ("good", "mua"), "untested condition"),
]

print(f"{'session':<34}{'units':>6}{'windows':>9}{'sig':>5}{'score':>8}{'p_session':>12}  verdict")
print("-" * 96)

for mouse, task, qual, note in CASES:
    sess = loaders.load_ebrains(mouse, task, quality=qual)
    _, z, _ = mo.build_activity(sess.spike_times, sess.t0, sess.t1)
    r = mo.windowed_sequence_test(z[::DECIM], BIN_S, n_surrogates=200)
    verdict = "SEQUENCES" if r["sequences"] else "none"
    print(f"{mouse + ' / ' + task:<34}{sess.n_units:>6}{r['n_windows']:>9}"
          f"{r['n_sig_windows']:>5}{r['frac_sig_windows']:>8.2f}{r['p_session']:>12.2e}"
          f"  {verdict:<10} <- {note}")

print("\nRequirement: BOTH wheel sessions must come out positive, or every null "
      "downstream is uninterpretable.")
