#!/usr/bin/env python3
"""
print_tables.py
---------------
Replicates Tables 1 and 2 from Soplata et al. (2023) by classifying
each simulation file's C-state percentage using a pre-trained HMM.

Workflow:
  1. Train (or load cached) HMM on a reference simulation file.
  2. Glob all Table 1 / Table 2 simulation files.
  3. Score each file in parallel using the frozen HMM.
  4. Aggregate scores by (gTC, gPY) parameter combination.
  5. Print formatted tables with mean ± std across seeds.

File naming convention assumed:
  sim_table1_gTC_0.0020_gPY_0.0020_seed_1.pkl
  sim_table2_gTC_0.0020_gPY_0.0030_seed_1.pkl
"""

import os
import re
import pickle
import argparse
import joblib
import numpy as np
from glob import glob
from collections import defaultdict
from scipy.signal import welch
from hmmlearn.hmm import GaussianHMM
from sklearn.preprocessing import StandardScaler


# ============================================================
# CONFIGURATION — edit paths and parameters here
# ============================================================

# Feature extraction parameters (physics-determined, not tunable)
# BIN_MS=5 → fs=200 Hz, Nyquist=100 Hz, fully resolves 8-20 Hz beta band
WINDOW_MS = 2000.0   # PSD window length; 2 s gives 0.5 Hz freq resolution
STEP_MS   = 100.0    # Stride between windows; controls temporal resolution
BIN_MS    = 5.0      # Spike histogram bin width in ms

# Burn-in period to discard from the start of every simulation (ms)
# Network needs time to settle into its propofol dynamics before scoring
BURN_IN   = 5000.0


# ============================================================
# STEP 1: FEATURE EXTRACTION
# ============================================================

def load_spikes(pkl_path, burn_in=BURN_IN):
    """
    Load spike data from a simulation pickle file, remove burn-in,
    and return post-burn-in spikes with t=0 at the end of burn-in.

    Parameters
    ----------
    pkl_path : str
        Path to the simulation .pkl file.
    burn_in : float
        Duration (ms) to discard from the start of the simulation.

    Returns
    -------
    tVec : np.ndarray
        Spike times in ms, shifted so t=0 is end of burn-in.
    idVec : np.ndarray
        Cell IDs corresponding to each spike.
    params : dict
        Network parameters (used for population sizes).
    """
    with open(pkl_path, 'rb') as f:
        data = pickle.load(f)

    spikes = data['spikes']
    params = data['network_params']

    # Shift time axis: t=0 is now the end of the burn-in period
    tVec  = np.array(spikes['tVec']) - burn_in
    idVec = np.array(spikes['idVec'])

    # Discard all spikes that occurred during burn-in
    valid = tVec >= 0
    return tVec[valid], idVec[valid], params


def build_ratio_series(pkl_path, window_ms=WINDOW_MS,
                       step_ms=STEP_MS, bin_ms=BIN_MS):
    """
    Convert a simulation file into a timeseries of delta/beta spectral
    power ratios using a sliding window Welch PSD.

    Each window's ratio = integral(PSD, 0.5-4 Hz) / integral(PSD, 8-20 Hz).
    High ratio → strong delta, synchronized (I-state candidate).
    Low ratio  → flat spectrum, desynchronized (C-state candidate).

    Parameters
    ----------
    pkl_path : str
        Path to the simulation .pkl file.
    window_ms, step_ms, bin_ms : float
        Window length, stride, and histogram bin width in ms.

    Returns
    -------
    ratios : np.ndarray, shape (n_windows,)
        Delta/beta power ratio for each window.
    times : np.ndarray, shape (n_windows,)
        Center time (ms, post-burn-in) of each window.
    """
    tVec, idVec, params = load_spikes(pkl_path)

    # Isolate cortical pyramidal cells (PYso + PYdr) — the classification signal
    py_end = params['Npyso'] + params['Npydr']
    py_t   = tVec[idVec < py_end]

    fs       = 1000.0 / bin_ms   # effective sampling rate of spike rate signal
    max_time = 30000.0  # We analyze up to 30s post-burn-in, matching the original paper's approach

    # Bin spikes into a continuous firing rate signal
    bins      = np.arange(0, max_time + bin_ms, bin_ms)
    counts, _ = np.histogram(py_t, bins=bins)

    w_bins = int(window_ms / bin_ms)   # window length in bins
    s_bins = int(step_ms   / bin_ms)   # stride in bins

    ratios, times = [], []

    for start in range(0, len(counts) - w_bins, s_bins):
        epoch = counts[start : start + w_bins]

        # Welch PSD: nperseg = full window → maximum frequency resolution
        freqs, psd = welch(epoch, fs=fs, nperseg=len(epoch))

        # Integrate power in each band using the trapezoidal rule
        lm    = (freqs >= 0.5) & (freqs <= 4.0)    # delta band
        hm    = (freqs >= 8.0) & (freqs <= 20.0)   # alpha/low-beta band
        low_p  = np.trapezoid(psd[lm], freqs[lm])
        high_p = np.trapezoid(psd[hm], freqs[hm])

        # Guard against division by zero in near-silent windows
        ratio = low_p / high_p if high_p > 0 else 100.0
        ratios.append(ratio)
        times.append(start * bin_ms + window_ms / 2)   # center of window

    return np.array(ratios), np.array(times)


# ============================================================
# STEP 2: TRAIN HMM (once) OR LOAD FROM CACHE
# ============================================================

def train_hmm(ref_path, cache_path, n_iter=50):
    """
    Fit a 2-state Gaussian HMM on the ratio timeseries of the reference
    simulation, then cache the result.  On subsequent calls the cache is
    loaded directly so the same decision boundary is applied to all files.

    The HMM learns:
      - means / covariances of the ratio distribution for each state
      - transition matrix A (dwell-time statistics)

    Polarity is resolved by convention: the state with the LOWER mean
    ratio is labelled C-state (desynchronized cortex → flat spectrum →
    low delta/beta ratio).

    Parameters
    ----------
    ref_path : str
        Path to the reference .pkl file used for training.
    cache_path : str
        Where to save/load the fitted model bundle.
    n_iter : int
        Maximum Baum-Welch iterations.

    Returns
    -------
    bundle : dict
        Keys: 'model', 'scaler', 'c_state_label', 'means'
    """
    if os.path.exists(cache_path):
        print(f"[HMM] Loading cached model from '{cache_path}'")
        return joblib.load(cache_path)

    print(f"[HMM] Training on reference file: '{ref_path}'")
    ratios, _ = build_ratio_series(ref_path)

    # Reshape to (n_samples, n_features=1) as required by hmmlearn
    X      = ratios.reshape(-1, 1)

    # StandardScaler improves Baum-Welch convergence when ratio range is large
    scaler = StandardScaler().fit(X)
    X_sc   = scaler.transform(X)

    model = GaussianHMM(
        n_components=2,
        covariance_type="full",
        n_iter=n_iter,
        random_state=42,
        tol=1e-4
    )
    model.fit(X_sc)

    # Inverse-transform learned means back to original ratio units
    means         = scaler.inverse_transform(model.means_).flatten()
    c_state_label = int(np.argmin(means))   # lower ratio = C-state

    # Log learned parameters for transparency
    print(f"[HMM] Learned state means (original units): {means.round(4)}")
    print(f"[HMM] C-state → label {c_state_label}  "
          f"(μ={means[c_state_label]:.4f})")
    print(f"[HMM] I-state → label {1-c_state_label}  "
          f"(μ={means[1-c_state_label]:.4f})")
    print(f"[HMM] Transition matrix:\n{model.transmat_.round(4)}")

    # Average dwell time per state: E[dwell] = 1 / (1 - A_ii) steps
    for i, lab in [(c_state_label, 'C'), (1 - c_state_label, 'I')]:
        dwell_steps = 1.0 / (1.0 - model.transmat_[i, i])
        dwell_s     = dwell_steps * STEP_MS / 1000.0
        print(f"[HMM] Expected {lab}-state dwell: "
              f"{dwell_steps:.1f} steps ≈ {dwell_s:.1f} s")

    bundle = {
        "model":         model,
        "scaler":        scaler,
        "c_state_label": c_state_label,
        "means":         means,
    }
    joblib.dump(bundle, cache_path)
    print(f"[HMM] Model cached to '{cache_path}'\n")
    return bundle


# ============================================================
# STEP 3: SCORE A SINGLE FILE
# ============================================================

def score_file(pkl_path, bundle):
    """
    Apply the pre-trained (frozen) HMM to one simulation file and
    return the percentage of time classified as C-state.

    The scaler and model are used read-only — no retraining occurs.
    This function is designed to be called in parallel via joblib.

    Parameters
    ----------
    pkl_path : str
        Path to a simulation .pkl file.
    bundle : dict
        Trained HMM bundle from train_hmm().

    Returns
    -------
    float
        Percentage of windows classified as C-state (0–100),
        or np.nan if the file contains no usable data.
    """
    try:
        ratios, _ = build_ratio_series(pkl_path)
        if len(ratios) == 0:
            return np.nan

        X_sc  = bundle["scaler"].transform(ratios.reshape(-1, 1))
        seq   = bundle["model"].predict(X_sc)   # Viterbi decoding
        return (seq == bundle["c_state_label"]).mean() * 100.0

    except Exception as e:
        print(f"  [WARN] Failed to score '{pkl_path}': {e}")
        return np.nan


# ============================================================
# STEP 4: PARSE PARAMETER VALUES FROM FILENAME
# ============================================================

def parse_filename(fname):
    """
    Extract gTC, gPY, and seed from a simulation filename.

    Expected pattern (anywhere in the filename):
        gTC_<float>_gPY_<float>_seed_<int>

    Example:
        sim_table1_gTC_0.0020_gPY_0.0030_seed_42.pkl
        → {'gTC': 0.002, 'gPY': 0.003, 'seed': 42}

    Returns None if the pattern is not found.
    """
    m = re.search(
        r'gTC_([\d.eE+\-]+)_gPY_([\d.eE+\-]+)_seed_(\d+)',
        fname
    )
    if m:
        return {
            'gTC':  float(m.group(1)),
            'gPY':  float(m.group(2)),
            'seed': int(m.group(3)),
        }
    return None


# ============================================================
# STEP 5: GLOB, SCORE IN PARALLEL, GROUP BY PARAMETERS
# ============================================================

def build_table(glob_pattern, bundle, n_jobs=-1):
    """
    Find all simulation files matching glob_pattern, score them in
    parallel, and group C-state percentages by (gTC, gPY) key.

    Parameters
    ----------
    glob_pattern : str
        Shell-style glob, e.g. 'results/data/tables_propofol/sim_table1_*.pkl'
    bundle : dict
        Trained HMM bundle from train_hmm().
    n_jobs : int
        Number of parallel workers (-1 = all CPU cores).

    Returns
    -------
    groups : dict
        Mapping (gTC, gPY) → list of C-state percentages (one per seed).
    skipped : int
        Number of files that could not be parsed or scored.
    """
    files = sorted(glob(glob_pattern))
    if not files:
        raise FileNotFoundError(
            f"No files found matching: '{glob_pattern}'\n"
            f"Check TABLE_DIR and filename prefix."
        )
    print(f"[Score] Found {len(files)} files. Scoring on {n_jobs} workers ...")

    # Parallel scoring — each call is independent, no shared state
    raw_scores = joblib.Parallel(n_jobs=n_jobs, verbose=2)(
        joblib.delayed(score_file)(f, bundle) for f in files
    )

    groups  = defaultdict(list)
    skipped = 0

    for fpath, score in zip(files, raw_scores):
        info = parse_filename(os.path.basename(fpath))

        # Skip files with unparseable names or failed scoring
        if info is None:
            print(f"  [WARN] Cannot parse filename: '{os.path.basename(fpath)}'")
            skipped += 1
            continue
        if np.isnan(score):
            print(f"  [WARN] NaN score for: '{os.path.basename(fpath)}'")
            skipped += 1
            continue

        groups[(info['gTC'], info['gPY'])].append(score)

    print(f"[Score] Done. {len(groups)} parameter combinations found. "
          f"{skipped} files skipped.\n")
    return groups, skipped


# ============================================================
# STEP 6: PRINT FORMATTED TABLES
# ============================================================

def print_table1(groups):
    """
    Print Table 1: each row is a (gTC, gPY) combination,
    columns show mean ± std C-state % and seed count.

    Table 1 in the paper varies a single parameter per row
    (70 independent seeds per row).
    """
    all_keys = sorted(groups.keys())
    if not all_keys:
        print("  [Table 1] No data to display.\n")
        return

    # Column widths
    w = {'gtc': 10, 'gpy': 10, 'mean': 10, 'std': 8, 'n': 5}
    sep = (f"  {'-'*w['gtc']}  {'-'*w['gpy']}  "
           f"{'-'*w['mean']}  {'-'*w['std']}  {'-'*w['n']}")

    print("=" * 60)
    print("  TABLE 1 — C-state % per parameter combination")
    print(f"  (mean ± std across seeds;  expected n=70 per row)")
    print("=" * 60)
    print(f"  {'gTC':>{w['gtc']}}  {'gPY':>{w['gpy']}}  "
          f"{'mean %':>{w['mean']}}  {'std %':>{w['std']}}  {'n':>{w['n']}}")
    print(sep)

    for key in all_keys:
        vals = groups[key]
        print(f"  {key[0]:>{w['gtc']}.4f}  {key[1]:>{w['gpy']}.4f}  "
              f"{np.mean(vals):>{w['mean']}.2f}  "
              f"{np.std(vals):>{w['std']}.2f}  "
              f"{len(vals):>{w['n']}}")

    print()


def print_table2(groups):
    """
    Print Table 2: a 2D grid with gTC as rows and gPY as columns.
    Each cell shows mean ± std C-state % across seeds.

    Table 2 in the paper sweeps both gTC and gPY
    (40 independent seeds per cell).
    """
    all_keys = sorted(groups.keys())
    if not all_keys:
        print("  [Table 2] No data to display.\n")
        return

    row_vals = sorted(set(k[0] for k in all_keys))   # unique gTC values
    col_vals = sorted(set(k[1] for k in all_keys))   # unique gPY values

    # Each cell: "XX.X±XX.X" — 10 chars wide
    cell_w   = 12
    label_w  = 10

    # Build header
    header = f"  {'gTC \\ gPY':>{label_w}}"
    for cv in col_vals:
        header += f"  {cv:>{cell_w}.4f}"

    divider = "  " + "-" * (label_w + (cell_w + 2) * len(col_vals))

    print("=" * len(divider))
    print("  TABLE 2 — C-state % grid: gTC (rows) × gPY (cols)")
    print(f"  (mean ± std across seeds;  expected n=40 per cell)")
    print("=" * len(divider))
    print(header)
    print(divider)

    for rv in row_vals:
        row_str = f"  {rv:>{label_w}.4f}"
        for cv in col_vals:
            vals = groups.get((rv, cv), [])
            if vals:
                # Format as "XX.X±XX.X" right-aligned within cell_w
                cell = f"{np.mean(vals):.1f}±{np.std(vals):.1f}"
            else:
                cell = "N/A"
            row_str += f"  {cell:>{cell_w}}"
        print(row_str)

    print()


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Generate Tables 1 and 2 from parameter sweeps.")
    parser.add_argument('--engine', type=str, default='NEURON-implementation',
                        choices=['NEURON-implementation', 'Brian2-implementation'],
                        help="Which simulator's data to analyze")
    parser.add_argument('--ref', type=str, default="sim_data_propofol_20260314_014344.pkl",
                        help="Filename of the reference simulation")
    args = parser.parse_args()

    # Dynamic paths based on engine choice
    ref_file_path = os.path.join(args.engine, "results", "data", "propofol", args.ref)
    table_dir_path = os.path.join(args.engine, "results", "data", "tables_propofol")
    model_cache = os.path.join(args.engine, "results", "models", "hmm_model.pkl")

    # Ensure the models directory exists for the cache
    os.makedirs(os.path.dirname(model_cache), exist_ok=True)

    # ----------------------------------------------------------
    # 1. Train HMM on reference file (or load from cache)
    #    The trained scaler + model are frozen after this step.
    #    All subsequent files are scored using the same boundary.
    # ----------------------------------------------------------
    bundle = train_hmm(ref_file_path, cache_path=model_cache)

    # ----------------------------------------------------------
    # 2. Table 1 — single-axis parameter sweep
    #    70 seeds per (gTC, gPY) row
    # ----------------------------------------------------------
    t1_pattern = os.path.join(table_dir_path, "sim_table1_*.pkl")
    t1_groups, t1_skipped = build_table(t1_pattern, bundle)
    print_table1(t1_groups)

    # ----------------------------------------------------------
    # 3. Table 2 — 2D parameter grid
    #    40 seeds per (gTC, gPY) cell
    # ----------------------------------------------------------
    t2_pattern = os.path.join(table_dir_path, "sim_table2_*.pkl")
    t2_groups, t2_skipped = build_table(t2_pattern, bundle)
    print_table2(t2_groups)

    # ----------------------------------------------------------
    # 4. Summary
    # ----------------------------------------------------------
    total_skipped = t1_skipped + t2_skipped
    if total_skipped > 0:
        print(f"[Summary] {total_skipped} files were skipped in total. "
              f"Check warnings above.")
    else:
        print("[Summary] All files scored successfully.")
