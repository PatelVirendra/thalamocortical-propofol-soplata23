import os
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt
import mne
from mne.time_frequency import tfr_array_multitaper
from scipy.signal import butter, filtfilt, spectrogram, decimate

# ==========================================
# 1. DATA LOADING & SETUP
# ==========================================
def load_latest_data(condition='propofol', base_dir=os.path.join('results', 'data')):
    """Finds and loads the most recent pickle file for a given condition."""
    target_dir = os.path.join(base_dir, condition)
    if not os.path.exists(target_dir):
        raise FileNotFoundError(f"Directory {target_dir} does not exist.")
    
    files = [f for f in os.listdir(target_dir) if f.endswith('.pkl')]
    if not files:
        raise FileNotFoundError(f"No pickle files found in {target_dir}.")
    
    # Get the latest file based on creation time
    latest_file = max([os.path.join(target_dir, f) for f in files], key=os.path.getctime)
    print(f"Loading data from: {latest_file}")
    
    with open(latest_file, 'rb') as f:
        data = pickle.load(f)
    return data

# ==========================================
# 2. SIGNAL PROCESSING (FILTERING)
# ==========================================
def bandpass_filter(data, lowcut, highcut, fs, order=2):
    """Applies a Butterworth bandpass filter (Order 2 as per paper)."""
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    y = filtfilt(b, a, data)
    return y

# ==========================================
# 3. PLOTTING FUNCTIONS
# ==========================================
def plot_zoomed_raster(spikes, params, save_dir, t_start, t_stop):
    """Replicates Fig 2C (right): Zoomed-in raster plot of thalamic alpha/low-beta."""
    tVec = np.array(spikes['tVec'])
    idVec = np.array(spikes['idVec'])
    
    # Filter spikes within the time window
    mask = (tVec >= t_start) & (tVec <= t_stop)
    t_zoom = tVec[mask]
    id_zoom = idVec[mask]
    
    plt.figure(figsize=(9.6, 10.8))
    plt.scatter(t_zoom, id_zoom, marker='|', s=10, color='black', alpha=0.6)
    
    # Add boundaries for cell types
    Npyso = params['Npyso']
    Npydr = params['Npydr']
    Ninh = params['Ninh']
    Ntc = params['Ntc']
    
    boundaries = [Npyso, Npyso+Npydr, Npyso+Npydr+Ninh, Npyso+Npydr+Ninh+Ntc]
    for b in boundaries:
        plt.axhline(b, color='gray', linestyle='--', linewidth=0.5)
        
    plt.xlim(t_start, t_stop)
    plt.xlabel('Time (ms)')
    plt.ylabel('Neuron ID')
    # plt.title('Zoomed-in Spike Rastergram (Propofol)')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'fig2c_zoomed_raster.png'), dpi=300)
    plt.close()

def plot_raw_eeg_window(time, eeg_signal, state_name, save_dir, t_start, t_stop):
    """Plots the raw, unfiltered EEG for a specific zoomed-in state window."""
    # 1. Calculate sampling frequency
    dt_ms = time[1] - time[0]
    fs = 1000.0 / dt_ms
    
    # 2. Apply standard clinical EEG filter (0.1 to 50 Hz) to remove spiking noise
    # This mimics the physical low-pass filtering of the skull/scalp
    eeg_clinical = bandpass_filter(eeg_signal, 0.1, 50.0, fs)
    
    # 3. Extract the specific time window
    mask = (time >= t_start) & (time <= t_stop)
    t_win = time[mask] / 1000.0
    eeg_win = eeg_clinical[mask]
    
    # 4. Plot
    plt.figure(figsize=(7, 3.5))
    plt.plot(t_win, eeg_win, color='black', linewidth=1)
    plt.xlim(t_win[0], t_win[-1])
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude (mV)')
    # plt.title(f'Raw Simulated EEG: {state_name}')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'fig2_raw_eeg_{state_name.lower()}.png'), dpi=300)
    plt.close()

def plot_eeg_spectrogram_multitaper(time, eeg_signal, save_dir):
    """Replicates Fig 2F: Propofol EEG Power Spectrum using Multitaper."""
    dt_ms = time[1] - time[0]
    fs_original = 1000.0 / dt_ms
    
    # 1. Filter EEG between 0.1 and 50 Hz
    eeg_filtered = bandpass_filter(eeg_signal, 0.1, 50.0, fs_original)
    
    # 2. Downsample the data to save memory and compute time
    # Downsampling from 10,000 Hz to 200 Hz (a factor of 50) is safe since our max freq is 50 Hz
    downsample_factor = int(fs_original / 200.0) 
    fs_new = fs_original / downsample_factor
    eeg_downsampled = decimate(eeg_filtered, downsample_factor, zero_phase=True)
    time_downsampled = time[::downsample_factor] / 1000.0 # Convert to seconds
    
    # 3. Format data for MNE: shape must be (n_epochs, n_channels, n_times)
    data = eeg_downsampled[np.newaxis, np.newaxis, :]
    
    # 4. Multitaper Parameters
    freqs = np.arange(0.1, 50.5, 0.5)
    n_cycles = freqs * 0.5          # ~0.25 sec window (sharper temporal resolution)

    power = tfr_array_multitaper(
        data,
        sfreq=fs_new,
        freqs=freqs,
        n_cycles=n_cycles,
        time_bandwidth=2.0,           # fewer tapers → sharper frequency bands
        output='power',
        n_jobs=1
    )

    spec_data = power[0, 0, :, :]
    spec_db = 10 * np.log10(spec_data)
    print(f"spec_db min: {spec_db.min():.1f}, max: {spec_db.max():.1f}, mean: {spec_db.mean():.1f}")
    vmin = np.percentile(spec_db, 5)   # ~ignore bottom 5%
    vmax = np.percentile(spec_db, 95)  # ~ignore top 5%
    print(f"Using vmin={vmin:.1f}, vmax={vmax:.1f}")

    # 5. Plotting
    plt.figure(figsize=(12, 4))
    
    # Use pcolormesh for a smooth, interpolated heat map
    plt.pcolormesh(time_downsampled, freqs, spec_db,
               cmap='jet', shading='gouraud',
               vmin=vmin, vmax=vmax)
    plt.ylim(0, 50)
    plt.xlabel('Time (sec)')
    plt.ylabel('Frequency (Hz)')
    # plt.title('Propofol EEG Power Spectrum (Multitaper)')
    plt.colorbar(label='Intensity (dBmV)')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'fig2f_eeg_spectrogram_multitaper.png'), dpi=300)
    plt.close()

def plot_eeg_states(time, eeg_signal, state_name, save_dir, t_start, t_stop):
    """Replicates Fig 3G and 3I: Filtered EEG for I-state and C-state."""
    dt_ms = time[1] - time[0]
    fs = 1000.0 / dt_ms
    
    # Apply the specific filters from the paper
    slow_wave = bandpass_filter(eeg_signal, 0.1, 1.0, fs)
    alpha_wave = bandpass_filter(eeg_signal, 8.0, 12.0, fs)
    low_beta_wave = bandpass_filter(eeg_signal, 8.0, 20.0, fs)
    
    # Extract window
    mask = (time >= t_start) & (time <= t_stop)
    t_win = time[mask] / 1000.0 # Convert to seconds for plotting
    
    plt.figure(figsize=(9, 7))
    plt.plot(t_win, slow_wave[mask], color='orange', label='Slow (0.1-1Hz)', linewidth=2)
    plt.plot(t_win, low_beta_wave[mask], color='red', label='Alpha/Low Beta (8-20Hz)', alpha=0.7)
    plt.plot(t_win, alpha_wave[mask], color='blue', label='Alpha (8-12Hz)', alpha=0.7)
    
    plt.xlim(t_win[0], t_win[-1])
    plt.xlabel('Time (sec)')
    plt.ylabel('Amplitude (mV)')
    # plt.title(f'{state_name} Simulated EEG')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'fig3_{state_name.lower()}_eeg.png'), dpi=300)
    plt.close()

def plot_individual_voltages(voltages, time, params, state_name, save_dir, t_start, t_stop):
    """Replicates Fig 3A and 3C: Voltage traces of individual cell types."""
    # Select one representative GID from each group
    gids = {
        'PYdr': params['Npyso'] + 50,             # Middle of PYdr
        'PYso': 50,                               # Middle of PYso
        'IN': params['Npyso'] + params['Npydr'] + 10,  # Middle of IN
        'TC': params['Npyso'] + params['Npydr'] + params['Ninh'] + 10, # Middle of TC
        'TRN': params['Npyso'] + params['Npydr'] + params['Ninh'] + params['Ntc'] + 10 # Middle of TRN
    }
    
    mask = (time >= t_start) & (time <= t_stop)
    t_win = time[mask] / 1000.0 # Seconds
    
    fig, axes = plt.subplots(5, 1, figsize=(9.6, 10.8), sharex=True)
    
    for ax, (cell_type, gid) in zip(axes, gids.items()):
        if gid in voltages:
            v_trace = voltages[gid][mask]
            ax.plot(t_win, v_trace, color='black', linewidth=0.5)
        ax.set_ylabel(f'{cell_type}\n[mV]')
        ax.set_ylim(-100, 50)
    
    axes[-1].set_xlabel('Time [sec]')
    # fig.suptitle(f'{state_name} Voltage Traces', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'fig3_{state_name.lower()}_voltages.png'), dpi=300)
    plt.close()

def plot_windowed_raster(spikes, params, save_dir, filename, t_start, t_stop):
    """Plots a raster for a specific time window using the exact styling from network.py."""
    tVecAll = np.array(spikes['tVec'])
    idVecAll = np.array(spikes['idVec'])
    
    # Filter for the specific time window
    mask = (tVecAll >= t_start) & (tVecAll <= t_stop)
    tVec_win = tVecAll[mask]
    idVec_win = idVecAll[mask]
    
    Npyso = params['Npyso']
    Npydr = params['Npydr']
    Ninh = params['Ninh']
    Ntc = params['Ntc']
    Nre = params['Nre']
    N = Npyso + Npydr + Ninh + Ntc + Nre

    # Indices for each cell type
    pyso_idx = [i for i, val in enumerate(idVec_win) if val < Npyso]
    pydr_idx = [i for i, val in enumerate(idVec_win) if Npyso <= val < Npyso + Npydr]
    inh_idx  = [i for i, val in enumerate(idVec_win) if Npyso + Npydr <= val < Npyso + Npydr + Ninh]
    tc_idx   = [i for i, val in enumerate(idVec_win) if Npyso + Npydr + Ninh <= val < Npyso + Npydr + Ninh + Ntc]
    re_idx   = [i for i, val in enumerate(idVec_win) if Npyso + Npydr + Ninh + Ntc <= val < N]

    cell_types = [
        ([tVec_win[i] for i in pyso_idx], [idVec_win[i] for i in pyso_idx], 'black', 'PYso'),
        ([tVec_win[i] for i in pydr_idx], [idVec_win[i] for i in pydr_idx], 'black', 'PYdr'),
        ([tVec_win[i] for i in inh_idx],  [idVec_win[i] for i in inh_idx],  'black', 'IN'),
        ([tVec_win[i] for i in tc_idx],   [idVec_win[i] for i in tc_idx],   'black', 'TC'),
        ([tVec_win[i] for i in re_idx],   [idVec_win[i] for i in re_idx],   'black', 'RE')
    ]

    id_ranges = [
        (0, Npyso), 
        (Npyso, Npyso + Npydr), 
        (Npyso + Npydr, Npyso + Npydr + Ninh), 
        (Npyso + Npydr + Ninh, Npyso + Npydr + Ninh + Ntc), 
        (Npyso + Npydr + Ninh + Ntc, N)
    ]

    fig, axes = plt.subplots(5, 1, figsize=(9.6, 10.8), sharex=True, constrained_layout=True)

    for ax, (tVec, ids, color, label), (id_min, id_max) in zip(axes, cell_types, id_ranges):
        spike_dict = {}
        for t, cell_id in zip(tVec, ids):
            t_sec = t / 1000.0  # Convert to seconds
            if cell_id not in spike_dict:
                spike_dict[cell_id] = []
            spike_dict[cell_id].append(t_sec)

        spike_times = [spike_dict[cell_id] for cell_id in sorted(spike_dict.keys())]
        cell_ids = sorted(spike_dict.keys())
        
        if spike_times:
            ax.eventplot(spike_times, lineoffsets=cell_ids, colors=color, linelengths=0.9, linewidths=0.9)
        else:
            ax.text(0.5, 0.5, 'No spikes', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, color='gray')

        ax.set_ylim(id_min, id_max)
        ax.invert_yaxis()
        ax.set_ylabel(f'{label}')

    axes[-1].set_xlabel('Time (s)')
    axes[-1].set_xlim(t_start / 1000.0, t_stop / 1000.0)

    plt.savefig(os.path.join(save_dir, filename), dpi=300, facecolor='white', edgecolor='none')
    plt.close()

def plot_full_simulation_overview(time, eeg_signal, voltages, params, save_dir):
    """Plots the full 30-second filtered EEG and full 30-second individual voltage traces."""
    t_sec = time / 1000.0
    
    # 1. Plot Full Simulated Filtered Total EEG (Fig 3H style)
    dt_ms = time[1] - time[0]
    fs = 1000.0 / dt_ms
    
    # Apply the specific filters from the paper
    slow_wave = bandpass_filter(eeg_signal, 0.1, 1.0, fs)
    alpha_wave = bandpass_filter(eeg_signal, 8.0, 12.0, fs)
    low_beta_wave = bandpass_filter(eeg_signal, 8.0, 20.0, fs)
    
    # Extract window for the entire simulation
    t_win = time / 1000.0 # Convert to seconds for plotting
    
    plt.figure(figsize=(9, 7))
    plt.plot(t_win, slow_wave, color='orange', label='Slow (0.1-1Hz)', linewidth=2)
    plt.plot(t_win, low_beta_wave, color='red', label='Alpha/Low Beta (8-20Hz)', alpha=0.7)
    plt.plot(t_win, alpha_wave, color='blue', label='Alpha (8-12Hz)', alpha=0.7)
    
    plt.xlim(t_win[0], t_win[-1])
    plt.xlabel('Time (sec)')
    plt.ylabel('Amplitude (mV)')
    # plt.title(f'Total Simulated Filtered EEG')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'fig3_total_filtered_eeg.png'), dpi=300)
    plt.close()

    # 2. Plot Full Individual Voltage Traces
    gids = {
        'PYso': params['Npyso'] // 2,
        'PYdr': params['Npyso'] + (params['Npydr'] // 2),
        'IN': params['Npyso'] + params['Npydr'] + (params['Ninh'] // 2),
        'TC': params['Npyso'] + params['Npydr'] + params['Ninh'] + (params['Ntc'] // 2),
        'RE': params['Npyso'] + params['Npydr'] + params['Ninh'] + params['Ntc'] + (params['Nre'] // 2)
    }
    
    fig, axes = plt.subplots(5, 1, figsize=(9.6, 10.8), sharex=True)
    for ax, (cell_type, gid) in zip(axes, gids.items()):
        if gid in voltages:
            ax.plot(t_sec, voltages[gid], color='black', linewidth=0.3)
        ax.set_ylabel(f'{cell_type}\n[mV]')
        ax.set_ylim(-100, 50)
    
    axes[-1].set_xlim(t_sec[0], t_sec[-1])
    axes[-1].set_xlabel('Time (s)')
    # fig.suptitle('Single Voltage Traces', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'full_simulation_voltages.png'), dpi=300)
    plt.close()

    # 3. Plot the multitaper spectrogram for the full propofol simulation (Fig 2F style)
    plot_eeg_spectrogram_multitaper(time, eeg_signal, save_dir)

# ==========================================
# 4. EXECUTION
# ==========================================
if __name__ == '__main__':
    condition_to_analyze = 'propofol'
    # 1. Setup Argument Parser
    parser = argparse.ArgumentParser(description="Analyze simulation EEG and raster data.")
    parser.add_argument('--engine', type=str, default='NEURON-implementation', 
                        choices=['NEURON-implementation', 'Brian2-implementation'],
                        help="Which simulator's data to analyze")
    args = parser.parse_args()
    
    # 2. Construct dynamic base directory for input data and figures output
    base_data_dir = os.path.join(args.engine, 'results', 'data')
    figures_dir = os.path.join(args.engine, 'results', 'figures', condition_to_analyze)
    
    os.makedirs(figures_dir, exist_ok=True)
    print(f"Figures will be saved to: {figures_dir}")
    
    # 3. Pass the dynamic path to your data loading function
    # (Update your function call to use base_dir=base_data_dir)
    try:
        data = load_latest_data(condition=condition_to_analyze, base_dir=base_data_dir)
    except FileNotFoundError as e:
        print(e)
        exit(1)

    params = data['network_params']
    spikes = data['spikes']
    voltages = data['voltages']
    time = data['eeg']['time']
    eeg_signal = data['eeg']['signal']

    # 4. Remove the DC Offset (Mean-centering)
    eeg_signal = eeg_signal - np.mean(eeg_signal)
    
    # 5. Generate and save plots
    # Define exact windows in milliseconds
    zoom_start, zoom_stop = 25000, 27000     # Fig 2C Zoomed propofol window
    i_start, i_stop = 18000, 23000         # I-state window
    c_start, c_stop = 10000, 15000         # C-state window

    print("Generating Full Simulation Overviews (EEG, Spectrogram & Voltages)...")
    plot_full_simulation_overview(time, eeg_signal, voltages, params, figures_dir)
    
    print("Generating Zoomed Raster Plot (23s to 25s)...")
    plot_windowed_raster(spikes, params, figures_dir, 'fig2c_zoomed_raster.png', zoom_start, zoom_stop)
    
    print("Generating I-State Plots...")
    plot_windowed_raster(spikes, params, figures_dir, 'fig3_istate_raster.png', i_start, i_stop)
    plot_raw_eeg_window(time, eeg_signal, 'I-state', figures_dir, t_start=18000, t_stop=20000)
    plot_eeg_states(time, eeg_signal, 'I-state', figures_dir, i_start, i_stop) # From previous code
    plot_individual_voltages(voltages, time, params, 'I-state', figures_dir, i_start, i_stop) # From previous code
    
    print("Generating C-State Plots...")
    plot_windowed_raster(spikes, params, figures_dir, 'fig3_cstate_raster.png', c_start, c_stop)
    plot_raw_eeg_window(time, eeg_signal, 'C-state', figures_dir, t_start=25000, t_stop=27000)
    plot_eeg_states(time, eeg_signal, 'C-state', figures_dir, c_start, c_stop) # From previous code
    plot_individual_voltages(voltages, time, params, 'C-state', figures_dir, c_start, c_stop) # From previous code

    print("All figures generated successfully!")
