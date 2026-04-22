#!/usr/bin/env python3
"""
analyze_perturbations.py
------------------------
Processes stimulus perturbation experiments to replicate the statistical
bar charts in Figures 4 and 5 of Soplata et al. (2023).

The key question per file:
  "After a stimulus is applied to a network in state X, how much of the
   post-stimulus window ends up in the target state Y?"

This is answered by:
  1. Building a delta/beta ratio timeseries from PY spikes (same feature
     extraction pipeline as the tables script).
  2. Applying a pre-trained frozen GaussianHMM to decode C/I state labels
     for each window (same model as the tables script — no retraining).
  3. Restricting the decoded sequence to the post-stimulus analysis window
     (stim_start + 200ms transient buffer → end_time).
  4. Computing the fraction of that window classified as the target state.
  5. Bucketing into All/Most/Half/Some/None categories for plotting.
"""

import os
import glob
import pickle
import argparse
import joblib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter
from scipy.signal import welch
from hmmlearn.hmm import GaussianHMM
from sklearn.preprocessing import StandardScaler


# ============================================================
# CONFIGURATION
# ============================================================

# Feature extraction parameters — must match tables script exactly
# so the ratio timeseries is on the same scale the HMM was trained on
WINDOW_MS = 2000.0   # PSD window length (ms); 2s → 0.5 Hz freq resolution
STEP_MS   = 100.0    # Stride between windows (ms)
BIN_MS    = 5.0      # Spike bin width (ms); fs=200 Hz, Nyquist=100 Hz

# Burn-in to discard from start of every simulation (ms)
BURN_IN   = 5000.0

# How long after stimulus onset to skip before scoring (ms)
# Allows the immediate stimulus transient to settle before classification
TRANSIENT_BUFFER_MS = 1000.0

# Minimum post-stimulus analysis window required to score a file (ms)
# Files with less usable time than this are skipped
MIN_ANALYSIS_MS = 2500.0


# ============================================================
# STEP 1: FEATURE EXTRACTION
# (identical pipeline to tables script)
# ============================================================

def load_spikes(pkl_path, burn_in=BURN_IN):
    """
    Load spike data, remove burn-in, return post-burn-in spikes
    with t=0 at the end of the burn-in period.

    Parameters
    ----------
    pkl_path : str
        Path to simulation .pkl file.
    burn_in : float
        Duration (ms) to discard from simulation start.

    Returns
    -------
    tVec : np.ndarray
        Spike times (ms), shifted so t=0 = end of burn-in.
    idVec : np.ndarray
        Cell IDs for each spike.
    params : dict
        Network parameters dict from the pickle.
    """
    with open(pkl_path, 'rb') as f:
        data = pickle.load(f)

    spikes = data['spikes']
    params = data['network_params']

    # Shift so burn-in end = t=0, then discard everything before
    tVec  = np.array(spikes['tVec']) - burn_in
    idVec = np.array(spikes['idVec'])
    valid = tVec >= 0

    return tVec[valid], idVec[valid], params


def build_ratio_series(pkl_path, window_ms=WINDOW_MS,
                       step_ms=STEP_MS, bin_ms=BIN_MS):
    """
    Convert spike data into a sliding-window delta/beta spectral ratio
    timeseries over the entire post-burn-in simulation.

    Each window ratio = integral(PSD, 0.5–4 Hz) / integral(PSD, 8–20 Hz).
    High ratio → strong delta, synchronized → I-state.
    Low ratio  → flat spectrum, desynchronized → C-state.

    Parameters
    ----------
    pkl_path : str
        Path to simulation .pkl file.
    window_ms, step_ms, bin_ms : float
        Window length, stride, bin width (all in ms).

    Returns
    -------
    ratios : np.ndarray, shape (n_windows,)
        Delta/beta ratio per window.
    times_ms : np.ndarray, shape (n_windows,)
        Center time (ms, post-burn-in) of each window.
    params : dict
        Network parameters (passed through for population sizing).
    """
    tVec, idVec, params = load_spikes(pkl_path)

    # Use only cortical PY cells as the classification signal
    py_end = params['Npyso'] + params['Npydr']
    py_t   = tVec[idVec < py_end]

    fs       = 1000.0 / bin_ms
    max_time = 10000.0  # We analyze up to 10000ms post-burn-in, matching the original paper's approach

    bins      = np.arange(0, max_time + bin_ms, bin_ms)
    counts, _ = np.histogram(py_t, bins=bins)

    w_bins = int(window_ms / bin_ms)
    s_bins = int(step_ms   / bin_ms)

    ratios, times_ms = [], []

    for start in range(0, len(counts) - w_bins, s_bins):
        epoch      = counts[start : start + w_bins]
        freqs, psd = welch(epoch, fs=fs, nperseg=len(epoch))

        lm    = (freqs >= 0.5) & (freqs <= 4.0)    # delta band
        hm    = (freqs >= 8.0) & (freqs <= 20.0)   # alpha/low-beta band
        low_p  = np.trapezoid(psd[lm], freqs[lm])
        high_p = np.trapezoid(psd[hm], freqs[hm])

        ratio = low_p / high_p if high_p > 0 else 100.0
        ratios.append(ratio)
        times_ms.append(start * bin_ms + window_ms / 2)

    return np.array(ratios), np.array(times_ms), params


# ============================================================
# STEP 2: TRAIN / LOAD HMM
# (identical to tables script — same model, same cache file)
# ============================================================

def train_hmm(ref_path, cache_path, n_iter=50):
    """
    Load the cached HMM trained on the reference file, or train it
    fresh if the cache does not exist.

    Using the same model as the tables script is intentional: it
    ensures that 'C-state' and 'I-state' labels mean the same thing
    across all analyses.

    Parameters
    ----------
    ref_path : str
        Reference .pkl file (only used if cache is absent).
    cache_path : str
        Path to saved model bundle (.pkl via joblib).
    n_iter : int
        Max Baum-Welch iterations (only used if retraining).

    Returns
    -------
    bundle : dict
        Keys: 'model', 'scaler', 'c_state_label', 'means'
    """
    if os.path.exists(cache_path):
        print(f"[HMM] Loading cached model from '{cache_path}'")
        return joblib.load(cache_path)

    print(f"[HMM] Cache not found. Training on '{ref_path}' ...")
    ratios, _, _ = build_ratio_series(ref_path)
    X      = ratios.reshape(-1, 1)
    scaler = StandardScaler().fit(X)
    X_sc   = scaler.transform(X)

    model = GaussianHMM(n_components=2, covariance_type="full",
                        n_iter=n_iter, random_state=42, tol=1e-4)
    model.fit(X_sc)

    # Lower mean ratio = flatter spectrum = desynchronized = C-state
    means         = scaler.inverse_transform(model.means_).flatten()
    c_state_label = int(np.argmin(means))

    print(f"[HMM] C-state label={c_state_label}, means={means.round(4)}")

    bundle = {"model": model, "scaler": scaler,
              "c_state_label": c_state_label, "means": means}
    joblib.dump(bundle, cache_path)
    print(f"[HMM] Cached to '{cache_path}'\n")
    return bundle


# ============================================================
# STEP 3: POST-STIMULUS STATE SCORING
# ============================================================

def score_post_stimulus(pkl_path, bundle,
                        actual_end_time=15000.0,
                        target_state='I-state'):
    """
    Classify the post-stimulus network state for a single simulation
    and return a categorical label (All/Most/Half/Some/None/SKIP).

    Strategy
    --------
    1. Build the full-simulation ratio timeseries and decode it with
       the frozen HMM → a sequence of C/I labels per 500ms window.
    2. Find the stimulus start time (from params, shifted by burn-in).
    3. Select only windows whose centers fall within:
           [stim_start + TRANSIENT_BUFFER_MS,  end_time]
    4. Compute fraction of those windows labelled as target_state.
    5. Map that fraction to All/Most/Half/Some/None via fixed thresholds
       that match the paper's categorical scheme.

    The key difference from the old TC-gap heuristic: state labels come
    from the same HMM used for Table 1/2 scoring, so the definition of
    C-state and I-state is perfectly consistent across all analyses.

    Parameters
    ----------
    pkl_path : str
        Path to simulation .pkl file.
    bundle : dict
        Trained HMM bundle from train_hmm().
    actual_end_time : float
        Absolute simulation end time in ms (before burn-in removal).
        Typically 15000 ms for perturbation experiments.
    target_state : str
        'I-state' or 'C-state' — what we are asking "how much of this?"

    Returns
    -------
    str
        One of: 'All', 'Most', 'Half', 'Some', 'None', 'SKIP'
        'SKIP' means the file had insufficient post-stimulus data.
    """
    # --- Build ratio timeseries and decode with frozen HMM ---
    ratios, times_ms, params = build_ratio_series(pkl_path)

    if len(ratios) == 0:
        return 'SKIP'

    X_sc      = bundle["scaler"].transform(ratios.reshape(-1, 1))
    state_seq = bundle["model"].predict(X_sc)   # Viterbi decoding
    # Boolean mask: True where the window is in C-state
    is_c = state_seq == bundle["c_state_label"]

    # --- Compute post-burn-in analysis window boundaries (ms) ---
    actual_stim_start = params.get('stim_start', 7500.0)
    stim_start_ms     = actual_stim_start  - BURN_IN
    end_time_ms       = actual_end_time    - BURN_IN

    # Skip the immediate transient after stimulus onset
    analysis_start_ms = stim_start_ms + TRANSIENT_BUFFER_MS

    # Require a minimum amount of usable post-stimulus time
    if (end_time_ms - analysis_start_ms) < MIN_ANALYSIS_MS:
        return 'SKIP'

    # --- Select windows whose centers fall inside the analysis window ---
    post_mask = (times_ms >= analysis_start_ms) & (times_ms <= end_time_ms)

    if post_mask.sum() == 0:
        return 'SKIP'

    post_c    = is_c[post_mask]     # C-state flags for post-stimulus windows
    post_i    = ~post_c             # I-state flags

    # --- Fraction of post-stimulus windows in target state ---
    if target_state == 'I-state':
        target_fraction = post_i.mean()
    elif target_state == 'C-state':
        target_fraction = post_c.mean()
    else:
        raise ValueError(f"Unknown target_state: '{target_state}'")

    # --- Bucket into categorical labels (Matched to Fig 4/5 Captions) ---
    if   target_fraction >= 0.95: return 'All'
    elif target_fraction >= 0.55: return 'Most'   # 55% to 95%
    elif target_fraction >= 0.45: return 'Half'   # ~50% (Tight band)
    elif target_fraction >= 0.10: return 'Some'   # 10% to 45%
    else:                         return 'None'   # <10%


# ============================================================
# STEP 4: DATA AGGREGATION
# ============================================================

def process_experiment_directory(condition, stim_type, target_state,
                                 bundle,
                                 base_dir=os.path.join('results', 'data'),
                                 actual_end_time=15000.0):
    """
    Score all simulation files in a condition directory that match
    the given stim_type, and aggregate results into a Counter.

    Parameters
    ----------
    condition : str
        Subdirectory name, e.g. 'propofol_mostly_c'.
    stim_type : str
        Value of params['stim_type'] to filter on, e.g. 'sync'/'desync'.
    target_state : str
        'I-state' or 'C-state'.
    bundle : dict
        Trained HMM bundle.
    base_dir : str
        Root directory containing condition subdirectories.
    actual_end_time : float
        Raw simulation end time in ms (before burn-in removal).

    Returns
    -------
    Counter
        Counts of each categorical label (All/Most/Half/Some/None).
    """
    search_path = os.path.join(base_dir, condition, "*.pkl")
    files       = glob.glob(search_path)

    if not files:
        print(f"  [WARN] No files found in: '{search_path}'")
        return Counter()

    print(f"  Scoring {len(files)} files in '{condition}' "
          f"(stim_type='{stim_type}', target='{target_state}') ...")

    scores  = []
    skipped = 0

    for fpath in sorted(files):
        with open(fpath, 'rb') as f:
            data = pickle.load(f)

        # Only process files with the requested stimulus type
        file_stim_type = data['network_params'].get('stim_type', 'none')
        if file_stim_type != stim_type:
            continue

        label = score_post_stimulus(
            fpath, bundle,
            actual_end_time=actual_end_time,
            target_state=target_state
        )

        if label == 'SKIP':
            skipped += 1
        else:
            scores.append(label)

    print(f"  -> {len(scores)} scored, {skipped} skipped.")
    return Counter(scores)


# ============================================================
# STEP 5: PLOTTING (unchanged from original)
# ============================================================

def plot_combined_stacked_bars(counts_sync, counts_desync,
                               title_sync, title_desync,
                               filename,
                               save_dir=os.path.join('results', 'figures')):
    """
    Produce a 1×2 stacked bar chart comparing sync vs desync stimulus
    responses, matching the layout of Figures 4D/E and 5D/E in the paper.

    Each bar shows the proportion of simulations falling into each
    categorical response label (All/Most/Half/Some/None).

    Parameters
    ----------
    counts_sync, counts_desync : Counter
        Category counts for synchronous and desynchronous stimuli.
    title_sync, title_desync : str
        Subplot titles (e.g. "D. Sync Input", "E. Desync Input").
    filename : str
        Output filename (saved inside save_dir).
    save_dir : str
        Directory to save the figure.
    """
    os.makedirs(save_dir, exist_ok=True)

    categories = ['All', 'Most', 'Half', 'Some', 'None']
    colors     = ['#1f77b4', '#d62728', '#ff7f0e', '#9467bd', '#2ca02c']

    def get_percentages(counts):
        total = sum(counts.values())
        if total == 0:
            return {cat: 0 for cat in categories}
        return {cat: (counts.get(cat, 0) / total) * 100
                for cat in categories}

    pct_sync   = get_percentages(counts_sync)
    pct_desync = get_percentages(counts_desync)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 5))
    plt.subplots_adjust(wspace=1.2)

    # --- Left bar: sync stimulus ---
    bottom = 0
    for cat, color in zip(categories, colors):
        val = pct_sync[cat]
        ax1.bar(0, val, bottom=bottom, color=color,
                edgecolor='white', width=0.5)
        bottom += val

    ax1.set_ylim(0, 100)
    ax1.set_xlim(-0.5, 0.5)
    ax1.set_xticks([])
    ax1.set_yticks([0, 25, 50, 75, 100])
    ax1.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
    ax1.set_ylabel('Proportion')
    ax1.set_title(title_sync, pad=15, fontsize=11, loc='center')

    # --- Right bar: desync stimulus ---
    bottom = 0
    for cat, color in zip(categories, colors):
        val = pct_desync[cat]
        ax2.bar(0, val, bottom=bottom, color=color,
                edgecolor='white', width=0.5)
        bottom += val

    ax2.set_ylim(0, 100)
    ax2.set_xlim(-0.5, 0.5)
    ax2.set_xticks([])
    ax2.set_yticks([0, 25, 50, 75, 100])
    ax2.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
    ax2.set_title(title_desync, pad=15, fontsize=11, loc='center')

    # --- Shared legend centred between the two bars ---
    patches = [mpatches.Patch(color=colors[i], label=categories[i])
               for i in range(len(categories))]
    fig.legend(handles=patches, title="Legend",
               loc='center', bbox_to_anchor=(0.5, 0.5))

    out_path = os.path.join(save_dir, filename)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  -> Saved: '{out_path}'")


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Analyze perturbation stimuli (Fig 4 & 5).")
    parser.add_argument('--engine', type=str, default='NEURON-implementation',
                        choices=['NEURON-implementation', 'Brian2-implementation'],
                        help="Which simulator's data to analyze")
    parser.add_argument('--ref', type=str, default="sim_data_propofol_20260314_014344.pkl",
                        help="Filename of the reference simulation")
    args = parser.parse_args()

    # Dynamic paths based on engine choice
    ref_file_path = os.path.join(args.engine, "results", "data", "propofol", args.ref)
    model_cache = os.path.join(args.engine, "results", "models", "hmm_model.pkl")
    base_data_dir = os.path.join(args.engine, "results", "data")
    figures_dir = os.path.join(args.engine, "results", "figures")

    # Ensure the models directory exists for the cache
    os.makedirs(os.path.dirname(model_cache), exist_ok=True)

    # ----------------------------------------------------------
    # 1. Load (or train) the shared HMM
    #    Same cache file as tables script → consistent state labels
    # ----------------------------------------------------------
    bundle = train_hmm(ref_file_path, cache_path=model_cache)

    # ----------------------------------------------------------
    # 2. Figure 4 — C-to-I transitions
    #    Network starts in mostly-C; stimulus tries to push it to I
    # ----------------------------------------------------------
    print("\n--- Figure 4: C-to-I transitions ---")

    fig4_sync = process_experiment_directory(
        condition='propofol_mostly_c', stim_type='sync',
        target_state='I-state', bundle=bundle, base_dir=base_data_dir
    )
    fig4_desync = process_experiment_directory(
        condition='propofol_mostly_c', stim_type='desync',
        target_state='I-state', bundle=bundle, base_dir=base_data_dir
    )

    if fig4_sync or fig4_desync:
        plot_combined_stacked_bars(
            fig4_sync, fig4_desync,
            title_sync="Sync Input", title_desync="Desync Input",
            filename="fig4_combined_bars.png", save_dir=figures_dir
        )
    else:
        print("  No data found for Figure 4.")

    # ----------------------------------------------------------
    # 3. Figure 5 — I-to-C transitions
    #    Network starts in mostly-I; stimulus tries to push it to C
    # ----------------------------------------------------------
    print("\n--- Figure 5: I-to-C transitions ---")

    fig5_sync = process_experiment_directory(
        condition='propofol_mostly_i', stim_type='sync',
        target_state='C-state', bundle=bundle, base_dir=base_data_dir
    )
    fig5_desync = process_experiment_directory(
        condition='propofol_mostly_i', stim_type='desync',
        target_state='C-state', bundle=bundle, base_dir=base_data_dir
    )

    if fig5_sync or fig5_desync:
        plot_combined_stacked_bars(
            fig5_sync, fig5_desync,
            title_sync="Sync Input", title_desync="Desync Input",
            filename="fig5_combined_bars.png", save_dir=figures_dir
        )
    else:
        print("  No data found for Figure 5.")

    print("\nAnalysis complete.")
