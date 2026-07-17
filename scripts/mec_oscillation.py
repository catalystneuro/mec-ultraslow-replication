"""Core library: detect minute-scale (ultraslow) oscillations and periodic
sequences in entorhinal populations.

Faithful to Gonzalo Cogno et al. (2024) Nature 625:338-344, Methods
(electrophysiology variant):

  1. bin spikes at 120 ms
  2. convolve each spike train with a Gaussian kernel, sigma = 5 s
  3. binarize at (mean + 1 SD) of each cell's smoothed trace
  4. autocorrelation of the binarized activity + power spectral density
  5. an "ultraslow oscillation" = a PSD peak below 0.1 Hz that is not
     present in shuffled versions of the same data
  6. PCA on the (time x cell) activity matrix; each cell's angle theta is the
     angle of its (PC1, PC2) loading vector; sorting cells by theta reveals
     periodic sequences

Deviations from the paper (documented deliberately):
  - PSD estimated with Welch's method rather than the FFT of the raw
    autocorrelation. Welch is a consistent (lower-variance) estimator of the
    same quantity; the raw periodogram has ~100% standard error at every
    frequency. The autocorrelogram is still computed and plotted.
  - Significance uses an inter-spike-interval (ISI) shuffle null, which is one
    of the two shuffles the paper used ("destroying the inter-event
    intervals"). It preserves each cell's spike count and ISI distribution but
    destroys slow rate modulation.
  - We additionally require the autocorrelogram to have a genuine secondary
    peak (a trough followed by a rebound). A monotonic slow drift -- from
    electrode drift, satiety, arousal decay -- also produces low-frequency PSD
    power but is NOT an oscillation. The paper's wheel-in-darkness recordings
    were far less exposed to this; datasets with task structure are, so we
    guard against it explicitly.
"""

import numpy as np
from scipy import signal, stats
from scipy.ndimage import gaussian_filter1d

BIN_SIZE = 0.12          # s, paper: 120 ms for Neuropixels
SIGMA = 5.0              # s, paper: Gaussian kernel width 5 s
BAND = (0.005, 0.1)      # Hz -> periods 10 s .. 200 s ("tens of seconds to minutes")


# --------------------------------------------------------------------------
# activity matrices
# --------------------------------------------------------------------------
def build_activity(spike_times, t0, t1, bin_size=BIN_SIZE, sigma=SIGMA):
    """Spike times -> (smoothed z-scored, binarized) activity matrices.

    Returns
    -------
    t     : (T,) bin centers, seconds
    z     : (T, N) Gaussian-smoothed, z-scored firing rate
    b     : (T, N) binarized at mean + 1 SD  (the paper's activity matrix)
    """
    edges = np.arange(t0, t1 + bin_size, bin_size)
    t = edges[:-1] + bin_size / 2
    n_cells = len(spike_times)
    z = np.zeros((len(t), n_cells))

    sigma_bins = sigma / bin_size
    for i, st in enumerate(spike_times):
        counts, _ = np.histogram(st, bins=edges)
        smoothed = gaussian_filter1d(counts.astype(float), sigma_bins, mode="reflect")
        sd = smoothed.std()
        z[:, i] = (smoothed - smoothed.mean()) / (sd if sd > 0 else 1.0)

    b = (z > 1.0).astype(float)   # z > 1 == raw > mean + 1 SD
    return t, z, b


# --------------------------------------------------------------------------
# spectral analysis
# --------------------------------------------------------------------------
def compute_psd(x, fs, nperseg_s=500.0):
    """Welch PSD of a single trace. nperseg chosen to resolve the ultraslow band."""
    nperseg = int(min(len(x), nperseg_s / (1 / fs)))
    f, p = signal.welch(x, fs=fs, nperseg=nperseg, noverlap=nperseg // 2,
                        detrend="constant", scaling="density")
    return f, p


def band_statistic(f, p, band=BAND):
    """Peak frequency in the ultraslow band and its share of total power.

    NOTE: kept for diagnostics only -- this statistic CANNOT detect the effect.
    The 5 s Gaussian kernel is a low-pass filter, so ~75% of the power sits in
    the ultraslow band for the real data AND for the ISI-shuffled surrogates
    alike (measured: 0.74 vs 0.77). Band power therefore carries almost no
    information about oscillatory structure. The real signature is a *narrowband
    peak* standing above the broadband 1/f-like background, which is what
    `oscillation_test` tests via the excess-power ratio below.
    """
    in_band = (f >= band[0]) & (f <= band[1])
    total = np.trapezoid(p, f)
    if total <= 0 or not in_band.any():
        return np.nan, 0.0
    idx = np.argmax(p[in_band])
    f_peak = f[in_band][idx]
    rel_power = p[in_band][idx] / total
    return f_peak, rel_power


def autocorr(x, max_lag_bins):
    """Normalized autocorrelation of a 1-D trace, lags 0..max_lag_bins."""
    x = x - x.mean()
    denom = np.dot(x, x)
    if denom == 0:
        return np.zeros(max_lag_bins + 1)
    full = signal.correlate(x, x, mode="full", method="fft")
    mid = len(x) - 1
    return full[mid:mid + max_lag_bins + 1] / denom


def has_secondary_peak(ac, bin_size=BIN_SIZE, min_period=10.0, max_period=200.0,
                       min_prominence=0.05):
    """True if the autocorrelogram dips below zero and then rebounds to a
    prominent peak -- the signature of an oscillation rather than a drift.
    Returns (bool, period_estimate_seconds)."""
    lo = int(min_period / bin_size)
    hi = min(int(max_period / bin_size), len(ac) - 1)
    if hi <= lo:
        return False, np.nan

    # require a trough below zero before the candidate peak
    peaks, props = signal.find_peaks(ac[:hi], prominence=min_prominence)
    peaks = peaks[peaks >= lo]
    if len(peaks) == 0:
        return False, np.nan
    first = peaks[0]
    if ac[:first].min() >= 0:      # never dipped below zero -> monotonic decay/drift
        return False, np.nan
    return True, first * bin_size


# --------------------------------------------------------------------------
# surrogates
# --------------------------------------------------------------------------
def isi_shuffle(st, t0, t1, rng):
    """Permute inter-spike intervals. Preserves spike count and ISI
    distribution; destroys slow rate modulation."""
    if len(st) < 3:
        return st.copy()
    isis = np.diff(np.sort(st))
    rng.shuffle(isis)
    out = st[0] + np.concatenate([[0.0], np.cumsum(isis)])
    return out[out <= t1]


def oscillation_test(spike_times, t0, t1, n_surrogates=200, band=BAND,
                     bin_size=BIN_SIZE, sigma=SIGMA, seed=0, progress=None):
    """Per-cell test for an ultraslow oscillation.

    Statistic: EXCESS POWER RATIO
        stat = max over f in band of  PSD_data(f) / mean_surrogate_PSD(f)

    Dividing by the surrogate mean spectrum whitens the broadband 1/f-like
    background that both data and surrogates share, so the statistic responds
    only to a narrowband peak rising above that background -- the actual
    signature of an oscillation. Taking the max over frequencies and comparing
    against the same max computed on each surrogate makes the test
    self-correcting for the search across frequencies.

    The null distribution is built by passing each surrogate through the
    identical statistic, using a leave-one-out surrogate mean so the observed
    and null statistics are computed the same way.

    Returns a dict of per-cell arrays:
        f_peak, excess_ratio, p_value, ac_period, has_peak
    """
    rng = np.random.default_rng(seed)
    fs = 1.0 / bin_size
    n_cells = len(spike_times)

    _, _, b_obs = build_activity(spike_times, t0, t1, bin_size, sigma)
    max_lag = int(300.0 / bin_size)

    f, _ = compute_psd(b_obs[:, 0], fs)
    fb = (f >= band[0]) & (f <= band[1])
    f_in = f[fb]

    psd_obs = np.zeros((n_cells, fb.sum()))
    ac_period = np.full(n_cells, np.nan)
    has_peak = np.zeros(n_cells, dtype=bool)
    for i in range(n_cells):
        _, p = compute_psd(b_obs[:, i], fs)
        psd_obs[i] = p[fb]
        ac = autocorr(b_obs[:, i], max_lag)
        has_peak[i], ac_period[i] = has_secondary_peak(ac, bin_size)

    # surrogate PSDs, in-band only (keeps memory small)
    psd_sur = np.zeros((n_surrogates, n_cells, fb.sum()))
    iterator = range(n_surrogates)
    if progress is not None:
        iterator = progress(iterator, desc="surrogates")
    for j in iterator:
        sur = [isi_shuffle(st, t0, t1, rng) for st in spike_times]
        _, _, b_sur = build_activity(sur, t0, t1, bin_size, sigma)
        for i in range(n_cells):
            _, p = compute_psd(b_sur[:, i], fs)
            psd_sur[j, i] = p[fb]

    eps = 1e-30
    tot = psd_sur.sum(axis=0)                       # (n_cells, n_f)
    mean_sur = tot / n_surrogates

    ratio_obs = psd_obs / np.maximum(mean_sur, eps)
    k_peak = np.argmax(ratio_obs, axis=1)
    excess_ratio = ratio_obs[np.arange(n_cells), k_peak]
    f_peak = f_in[k_peak]

    # null: same statistic on each surrogate, leave-one-out mean
    n_ge = np.zeros(n_cells)
    for j in range(n_surrogates):
        loo = (tot - psd_sur[j]) / (n_surrogates - 1)
        stat_j = (psd_sur[j] / np.maximum(loo, eps)).max(axis=1)
        n_ge += (stat_j >= excess_ratio)

    p_value = (1.0 + n_ge) / (1.0 + n_surrogates)
    return dict(f_peak=f_peak, excess_ratio=excess_ratio, p_value=p_value,
                ac_period=ac_period, has_peak=has_peak)


def benjamini_hochberg(p, alpha=0.05):
    """BH-FDR. Returns boolean array of rejections."""
    p = np.asarray(p)
    n = len(p)
    order = np.argsort(p)
    thresh = alpha * (np.arange(1, n + 1) / n)
    passed = p[order] <= thresh
    rej = np.zeros(n, dtype=bool)
    if passed.any():
        k = np.max(np.where(passed)[0])
        rej[order[:k + 1]] = True
    return rej


# --------------------------------------------------------------------------
# population sequences (paper Fig. 3)
# --------------------------------------------------------------------------
def pca_sequence(activity):
    """PCA on the (time x cell) activity matrix.

    Returns
    -------
    scores  : (T, 2) projection onto PC1, PC2  -- the oscillation trajectory
    theta   : (N,) each cell's angle from its (PC1, PC2) loading vector
    order   : (N,) cell indices sorted by descending theta (the sequence order)
    var_exp : (2,) fraction of variance explained by PC1, PC2
    phase   : (T,) instantaneous oscillation phase, atan2(PC2, PC1)
    """
    x = activity - activity.mean(axis=0, keepdims=True)
    # cells with zero variance carry no loading
    sd = x.std(axis=0)
    sd[sd == 0] = 1.0

    u, s, vt = np.linalg.svd(x, full_matrices=False)
    var_exp = (s ** 2) / np.sum(s ** 2)

    loadings = vt[:2].T                        # (N, 2)
    scores = u[:, :2] * s[:2]                  # (T, 2)
    theta = np.arctan2(loadings[:, 1], loadings[:, 0])
    order = np.argsort(-theta)
    phase = np.arctan2(scores[:, 1], scores[:, 0])
    return scores, theta, order, var_exp[:2], phase


def phase_locking(activity, phase, n_bins=36):
    """Per-cell phase tuning w.r.t. the population oscillation phase.

    Returns preferred phase, mean resultant length (Rayleigh R), and the
    Rayleigh p-value, weighting each time bin by that cell's activity.
    """
    n_cells = activity.shape[1]
    pref = np.full(n_cells, np.nan)
    mrl = np.zeros(n_cells)
    pval = np.ones(n_cells)

    for i in range(n_cells):
        w = activity[:, i]
        if w.sum() <= 0:
            continue
        # circular mean weighted by activity
        c = np.sum(w * np.cos(phase))
        s = np.sum(w * np.sin(phase))
        r = np.sqrt(c ** 2 + s ** 2) / w.sum()
        pref[i] = np.arctan2(s, c)
        mrl[i] = r
        # Rayleigh test with effective n = number of active bins
        n_eff = w.sum() ** 2 / np.sum(w ** 2)
        z = n_eff * r ** 2
        pval[i] = np.exp(-z) * (1 + (2 * z - z ** 2) / (4 * n_eff))
    return pref, mrl, pval


def acg_spectrum(counts, window=560.0, bin_size=BIN_SIZE, fs=None):
    """PSD of the z-scored spike-train autocorrelogram (detector 'C').

    Chosen over the PSD of the activity trace because the bake-off against the
    paper's own wheel data showed it is the only variant that recovers the
    effect (64% of cells vs 0-2%). The autocorrelogram averages over the whole
    session, which survives the oscillation's non-stationary phase; Welch on the
    raw trace averages segments and smears the peak away.
    """
    fs = fs or 1.0 / bin_size
    nlag = int(window / bin_size)
    x = counts - counts.mean()
    full = signal.correlate(x, x, mode="full", method="fft")
    mid = len(x) - 1
    a = full[max(0, mid - nlag): mid + nlag + 1]
    sd = a.std()
    a = (a - a.mean()) / (sd if sd > 0 else 1.0)
    return signal.welch(a, fs, window="hamming", nperseg=min(len(a), 8196),
                        detrend="constant")


def _bin_counts(spike_times, t0, t1, bin_size=BIN_SIZE):
    edges = np.arange(t0, t1 + bin_size, bin_size)
    return np.array([np.histogram(st, bins=edges)[0].astype(float)
                     for st in spike_times])


def single_cell_test(spike_times, t0, t1, n_surrogates=100, band=BAND,
                     bin_size=BIN_SIZE, seed=0, progress=None):
    """Per-cell ultraslow rhythmicity: excess autocorrelogram-spectrum power in
    the ultraslow band, against an ISI-shuffle null.

    The null is POOLED across cells: the excess-power ratio is normalized, hence
    exchangeable across cells under H0. Pooling lifts the p-value resolution from
    1/(n_surr+1) to ~1/(n_surr * n_cells), which is what makes BH-FDR across
    hundreds of cells achievable at all. Without it, the test returns 0
    significant cells on the paper's own data no matter how strong the effect.
    """
    rng = np.random.default_rng(seed)
    n_cells = len(spike_times)
    fs = 1.0 / bin_size

    cnt = _bin_counts(spike_times, t0, t1, bin_size)
    f, _ = acg_spectrum(cnt[0], fs=fs)
    fb = (f >= band[0]) & (f <= band[1])
    f_in = f[fb]

    psd_obs = np.array([acg_spectrum(c, fs=fs)[1][fb] for c in cnt])

    psd_sur = np.zeros((n_surrogates, n_cells, fb.sum()))
    it = range(n_surrogates)
    if progress is not None:
        it = progress(it, desc="ISI surrogates", leave=False)
    for j in it:
        sur = [isi_shuffle(st, t0, t1, rng) for st in spike_times]
        cs = _bin_counts(sur, t0, t1, bin_size)
        for i in range(n_cells):
            psd_sur[j, i] = acg_spectrum(cs[i], fs=fs)[1][fb]

    eps = 1e-30
    tot = psd_sur.sum(axis=0)
    mean_sur = tot / n_surrogates
    ratio = psd_obs / np.maximum(mean_sur, eps)
    k = np.argmax(ratio, axis=1)
    stat = ratio[np.arange(n_cells), k]
    f_peak = f_in[k]

    null = np.concatenate([
        (psd_sur[j] / np.maximum((tot - psd_sur[j]) / (n_surrogates - 1), eps)).max(axis=1)
        for j in range(n_surrogates)
    ])
    p = np.array([(1 + (null >= s).sum()) / (1 + len(null)) for s in stat])
    sig = benjamini_hochberg(p, 0.05)
    return dict(f_peak=f_peak, stat=stat, p_value=p, sig=sig, null=null)


# --------------------------------------------------------------------------
# population sequence test
# --------------------------------------------------------------------------
def _pop_metrics(z):
    """Two population statistics from the (time x cell) activity matrix.

    pc12_var : fraction of variance in the leading 2-D subspace. A single shared
               oscillation concentrates variance into one plane.
    rotation : consistency of phase advance of the PC1/PC2 trajectory,
               |mean(sign(dphi))|, in [0, 1]. A genuine oscillation sweeps the
               plane in a consistent direction; independently modulated cells
               produce a trajectory that wanders.
    """
    x = z - z.mean(axis=0, keepdims=True)
    c = np.cov(x, rowvar=False)
    w, v = np.linalg.eigh(c)
    idx = np.argsort(-w)
    w, v = w[idx], v[:, idx]
    pc12_var = float((w[0] + w[1]) / w.sum()) if w.sum() > 0 else 0.0

    scores = x @ v[:, :2]
    phase = np.arctan2(scores[:, 1], scores[:, 0])
    dphi = np.angle(np.exp(1j * np.diff(phase)))
    rotation = float(abs(np.mean(np.sign(dphi))))

    theta = np.arctan2(v[:, 1], v[:, 0])
    return pc12_var, rotation, theta, phase, scores, w[:2] / w.sum()


def population_sequence_test(z, n_surrogates=200, seed=0, progress=None):
    """Is there a SHARED periodic sequence, beyond independent slow modulation?

    Null = independent CIRCULAR SHIFT of each cell's activity trace. A circular
    shift leaves each cell's power spectrum and autocorrelation exactly intact --
    every cell keeps whatever ultraslow rhythmicity it had -- while destroying
    the timing relationships BETWEEN cells. So this null isolates precisely the
    paper's claim: not "do cells modulate slowly" (grid cells do that in an open
    field for trivial spatial reasons) but "do they ride a common oscillation,
    at staggered phases, forming a repeating sequence".
    """
    rng = np.random.default_rng(seed)
    T = z.shape[0]
    pc12, rot, theta, phase, scores, ve = _pop_metrics(z)

    null_pc12, null_rot = [], []
    it = range(n_surrogates)
    if progress is not None:
        it = progress(it, desc="circular shifts", leave=False)
    for _ in it:
        zs = np.column_stack([np.roll(z[:, i], int(rng.integers(T)))
                              for i in range(z.shape[1])])
        a, b, *_ = _pop_metrics(zs)
        null_pc12.append(a)
        null_rot.append(b)

    null_pc12 = np.array(null_pc12)
    null_rot = np.array(null_rot)
    p_pc12 = (1 + (null_pc12 >= pc12).sum()) / (1 + n_surrogates)
    p_rot = (1 + (null_rot >= rot).sum()) / (1 + n_surrogates)

    return dict(pc12_var=pc12, rotation=rot, theta=theta, phase=phase,
                scores=scores, var_exp=ve,
                p_pc12=float(p_pc12), p_rotation=float(p_rot),
                null_pc12_mean=float(null_pc12.mean()),
                null_rot_mean=float(null_rot.mean()),
                pc12_z=float((pc12 - null_pc12.mean()) / (null_pc12.std() + 1e-12)),
                rot_z=float((rot - null_rot.mean()) / (null_rot.std() + 1e-12)))


def windowed_sequence_test(z, bin_s, window_s=300.0, step_s=100.0,
                           n_surrogates=200, seed=0):
    """Sequence test that survives INTERMITTENCY.

    Why this exists: a whole-session rotation index averages the sequence epochs
    together with the quiet stretches, and the effect washes out. The paper is
    explicit that this happens -- "the periodic sequences were more salient in some
    subsets of the sessions than others" -- and it sorts on a 300 s subset for its
    own Neuropixels examples. The whole-session test duly FAILED on the paper's
    second Neuropixels mouse (102335, rot_z = -5.2), which is a false negative on
    known-positive data.

    So: slide a window, score each window against a circular-shift null, and ask
    whether SIGNIFICANT WINDOWS ARE MORE COMMON THAN CHANCE (binomial test). This
    is the same idea as the paper's 'oscillation score' (fraction of the session
    spent oscillating) but with a per-window null attached.

    Returns per-window statistics plus a session-level binomial p-value.
    """
    rng = np.random.default_rng(seed)
    T, N = z.shape
    w = int(window_s / bin_s)
    step = int(step_s / bin_s)
    if T < w:
        w, step = T, max(1, T // 2)

    starts = list(range(0, max(1, T - w + 1), step))
    rot_obs = np.zeros(len(starts))
    pc_obs = np.zeros(len(starts))
    for i, s in enumerate(starts):
        pc, rot, *_ = _pop_metrics(z[s:s + w])
        rot_obs[i], pc_obs[i] = rot, pc

    # null: circular-shift the FULL traces (preserving each cell's spectrum),
    # then score the same windows
    rot_null = np.zeros((n_surrogates, len(starts)))
    for j in range(n_surrogates):
        zs = np.column_stack([np.roll(z[:, i], int(rng.integers(T))) for i in range(N)])
        for i, s in enumerate(starts):
            _, r, *_ = _pop_metrics(zs[s:s + w])
            rot_null[j, i] = r

    # per-window p, then a session-level binomial test on the count of hits
    p_win = np.array([(1 + (rot_null[:, i] >= rot_obs[i]).sum()) / (1 + n_surrogates)
                      for i in range(len(starts))])
    n_sig = int((p_win < 0.05).sum())
    n_win = len(starts)
    p_session = float(stats.binomtest(n_sig, n_win, 0.05, alternative="greater").pvalue)

    # ---- the same test on INDEPENDENT (non-overlapping) windows.
    # The binomial above assumes the windows are independent, but with the defaults
    # (300 s window, 100 s step) adjacent windows share 2/3 of their samples, so it
    # counts ~3x more evidence than exists and is anti-conservative. Non-overlapping
    # windows are exactly every k-th window of this grid, so the corrected statistic
    # costs nothing extra to compute.
    #
    # This matters: on the paper's own mouse 102335 the overlapping test gives
    # 3/13 windows, p=0.025 (significant) while the independent test gives 1/5,
    # p=0.23 (not) -- and that session's median observed rotation is BELOW its median
    # null. The correction weakens POSITIVES (and the weaker positive control) while
    # leaving NULLS untouched or stronger, since a liberal test that found nothing
    # would find nothing under a stricter one.
    k = max(1, int(round(window_s / step_s)))
    sub = np.arange(0, len(starts), k)
    n_sig_i = int((p_win[sub] < 0.05).sum())
    n_win_i = len(sub)
    p_session_i = float(stats.binomtest(n_sig_i, n_win_i, 0.05,
                                        alternative="greater").pvalue)

    return dict(
        n_windows=n_win, n_sig_windows=n_sig,
        frac_sig_windows=float(n_sig / n_win),      # ~ the paper's 'oscillation score'
        p_session=p_session,
        n_windows_indep=n_win_i, n_sig_windows_indep=n_sig_i,
        frac_sig_windows_indep=float(n_sig_i / n_win_i),
        p_session_indep=p_session_i,
        best_rot_z=float((rot_obs - rot_null.mean(0)).max()
                         / (rot_null.std(0).mean() + 1e-12)),
        median_rot=float(np.median(rot_obs)),
        median_rot_null=float(np.median(rot_null)),
        sequences=bool(p_session < 0.05),
        sequences_indep=bool(p_session_i < 0.05),
    )


def sequence_score(pref_phase, mrl, sig_mask):
    """How uniformly do preferred phases tile the oscillation cycle?

    A *sequence* requires cells to fire at staggered phases spanning the whole
    cycle. A synchronous population burst would instead concentrate every cell
    at one phase. We therefore test the preferred phases of the significantly
    locked cells for UNIFORMITY: high p (non-significant Rayleigh on the
    phase distribution) = uniformly tiled = sequence.
    """
    ph = pref_phase[sig_mask & np.isfinite(pref_phase)]
    if len(ph) < 5:
        return dict(n=len(ph), uniformity_p=np.nan, spread=np.nan)
    c, s = np.cos(ph).sum(), np.sin(ph).sum()
    r = np.sqrt(c ** 2 + s ** 2) / len(ph)          # low r = uniform = sequence
    z = len(ph) * r ** 2
    p_rayleigh = np.exp(-z)                          # small p = clustered = NOT a sequence
    return dict(n=len(ph), uniformity_p=float(p_rayleigh), spread=float(1 - r))
