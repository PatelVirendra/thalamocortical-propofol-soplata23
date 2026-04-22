"""
network.py — Brian2 Thalamocortical Network
============================================
Defines the Net class, which assembles all NeuronGroups, Synapses, and
monitors; runs the simulation; and provides methods to plot raster plots,
print spike rates, and save data to disk.

Mirrors the structure of the NEURON network.py as closely as possible so
that post-processing and analysis scripts (HMM classifier, table scripts,
etc.) can operate on identical data formats from both simulators.
"""

import os
import pickle
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use('Agg')           # non-interactive backend for cluster use
import matplotlib.pyplot as plt
import brian2 as b2

import config
import cells
import synapses

class Net:
    """
    Full thalamocortical network for propofol anaesthesia modelling.

    Usage
    -----
    >>> net = Net()
    >>> net.run()
    >>> net.plotRaster()
    >>> net.printSpikeRates()
    >>> net.saveData()
    """

    def __init__(self):
        self.params = config.current_state_params
        self._build()

    # =========================================================================
    # NETWORK CONSTRUCTION
    # =========================================================================

    def _build(self):
        """Instantiate all groups, synapses, and monitors."""
        p = self.params

        # --- Cell groups ---
        (self.PYdr_group,
         self.poisson_stim_PYdr,
         self.poisson_syn_PYdr) = cells.make_PYdr_group(p)

        self.PYso_group = cells.make_PYso_group(p)

        self.IN_group = cells.make_IN_group(p)

        if p['include_thalamus']:
            (self.TC_group,
             self.poisson_stim_TC,
             self.poisson_syn_TC) = cells.make_TC_group(p)

            (self.RE_group,
             self.poisson_stim_RE,
             self.poisson_syn_RE) = cells.make_RE_group(p)

        # --- Compartmental / non-synaptic connections ---
        self.ICOM_PYdr_PYso, self.ICOM_PYso_PYdr = synapses.make_ICOM(
            self.PYdr_group, self.PYso_group)
        self.IKNa_PYso_PYdr = synapses.make_IKNa(
            self.PYdr_group, self.PYso_group)

        # --- Intracortical synapses ---
        self.IAMPA_PYso_PYdr = synapses.make_IAMPA_PYso_PYdr(
            self.PYso_group, self.PYdr_group, p)
        self.INMDA_PYso_PYdr = synapses.make_INMDA_PYso_PYdr(
            self.PYso_group, self.PYdr_group, p)
        self.IAMPA_PYso_IN   = synapses.make_IAMPA_PYso_IN(
            self.PYso_group, self.IN_group, p)
        self.INMDA_PYso_IN   = synapses.make_INMDA_PYso_IN(
            self.PYso_group, self.IN_group, p)
        self.IGABAA_IN_PYso  = synapses.make_IGABAA_IN_PYso(
            self.IN_group, self.PYso_group, p)
        self.IGABAA_IN_IN    = synapses.make_IGABAA_IN_IN(
            self.IN_group, p)

        # --- Intrathalamic & thalamocortical synapses ---
        if p['include_thalamus']:
            self.IAMPA_TC_RE    = synapses.make_IAMPA_TC_RE(
                self.TC_group, self.RE_group, p)
            self.IGABAA_RE_RE   = synapses.make_IGABAA_RE_RE(
                self.RE_group, p)
            self.IGABAA_RE_TC   = synapses.make_IGABAA_RE_TC(
                self.RE_group, self.TC_group, p)
            self.IGABAB_RE_TC   = synapses.make_IGABAB_RE_TC(
                self.RE_group, self.TC_group, p)

        if p['include_thalamocortical_connections']:
            self.IAMPA_TC_PYdr  = synapses.make_IAMPA_TC_PYdr(
                self.TC_group, self.PYdr_group, p)
            self.IAMPA_TC_IN    = synapses.make_IAMPA_TC_IN(
                self.TC_group, self.IN_group, p)
            self.IAMPA_PYso_TC  = synapses.make_IAMPA_PYso_TC(
                self.PYso_group, self.TC_group, p)
            self.IAMPA_PYso_RE  = synapses.make_IAMPA_PYso_RE(
                self.PYso_group, self.RE_group, p)

        # --- Spike and state monitors ---
        self._create_monitors()

        # --- Assemble the Brian2 Network ---
        self._assemble_network()

    def _create_monitors(self):
        """Create spike monitors and (optionally) voltage state monitors."""
        downsample_clock = b2.Clock(dt=config.record_interval)

        # Spike monitors — always created
        self.PYdr_spikemon = b2.SpikeMonitor(self.PYdr_group, name='PYdr_spikemon')
        self.PYso_spikemon = b2.SpikeMonitor(self.PYso_group, name='PYso_spikemon')
        self.IN_spikemon   = b2.SpikeMonitor(self.IN_group,   name='IN_spikemon')
        if self.params['include_thalamus']:
            self.TC_spikemon = b2.SpikeMonitor(self.TC_group, name='TC_spikemon')
            self.RE_spikemon = b2.SpikeMonitor(self.RE_group, name='RE_spikemon')

        # Voltage state monitors — created only when --record_voltages is set
        if config.record_voltages:
            self.PYdr_monitor = b2.StateMonitor(
                self.PYdr_group, 'v', record=True,
                clock=downsample_clock, name='PYdr_monitor')
            self.PYso_monitor = b2.StateMonitor(
                self.PYso_group, 'v', record=True,
                clock=downsample_clock, name='PYso_monitor')
            self.IN_monitor = b2.StateMonitor(
                self.IN_group, 'v', record=True,
                clock=downsample_clock, name='IN_monitor')
            if self.params['include_thalamus']:
                self.TC_monitor = b2.StateMonitor(
                    self.TC_group, 'v', record=True,
                    clock=downsample_clock, name='TC_monitor')
                self.RE_monitor = b2.StateMonitor(
                    self.RE_group, 'v', record=True,
                    clock=downsample_clock, name='RE_monitor')

    def _assemble_network(self):
        """Collect all Brian2 objects into a single Network instance."""
        objs = [
            # Cell groups
            self.PYdr_group, self.PYso_group, self.IN_group,
            # Poisson drives
            self.poisson_stim_PYdr, self.poisson_syn_PYdr,
            # Compartmental connections
            self.ICOM_PYdr_PYso, self.ICOM_PYso_PYdr, self.IKNa_PYso_PYdr,
            # Intracortical synapses
            self.IAMPA_PYso_PYdr, self.INMDA_PYso_PYdr,
            self.IAMPA_PYso_IN,   self.INMDA_PYso_IN,
            self.IGABAA_IN_PYso,  self.IGABAA_IN_IN,
            # Spike monitors
            self.PYdr_spikemon, self.PYso_spikemon, self.IN_spikemon,
        ]

        if self.params['include_thalamus']:
            objs += [
                self.TC_group, self.RE_group,
                self.poisson_stim_TC, self.poisson_syn_TC,
                self.poisson_stim_RE, self.poisson_syn_RE,
                self.IAMPA_TC_RE,   self.IGABAA_RE_RE,
                self.IGABAA_RE_TC,  self.IGABAB_RE_TC,
                self.TC_spikemon,   self.RE_spikemon,
            ]

        if self.params['include_thalamocortical_connections']:
            objs += [
                self.IAMPA_TC_PYdr, self.IAMPA_TC_IN,
                self.IAMPA_PYso_TC, self.IAMPA_PYso_RE,
            ]

        if config.record_voltages:
            objs += [self.PYdr_monitor, self.PYso_monitor, self.IN_monitor]
            if self.params['include_thalamus']:
                objs += [self.TC_monitor, self.RE_monitor]

        self.network = b2.Network(*objs)

    # =========================================================================
    # SIMULATION EXECUTION
    # =========================================================================

    def run(self):
        """
        Run the simulation via the C++ standalone device.
 
        b2.seed() is called just before network.run() so that Brian2's
        internal C++ PRNG — including the PoissonGroup spike trains — is
        seeded deterministically.  Combined with the numpy-based
        get_independent_draws() initialisation in cells.py / synapses.py,
        this guarantees bit-identical results across runs with the same
        --seed value.
        """
        print(f"\nRunning simulation for {config.duration / b2.ms:.0f} ms ...")
 
        # --- Seed Brian2's C++ RNG with a fixed constant ---
        # In NEURON, NetStim with noise=1 uses NEURON's internal Random[0],
        # which is always initialised to the same default state at startup
        # regardless of the --seed argument.  This means Poisson noise is
        # IDENTICAL across all independent trials in NEURON; only the cell/
        # synapse initial conditions (set via get_independent_draws) vary.
        # To match this behaviour in Brian2, we seed the C++ PRNG with a
        # fixed constant (0) rather than config.randSeed.  Using randSeed
        # here would incorrectly make the Poisson spike trains differ between
        # trials, which NEURON does not do. However, if it is somehow desired
        # to have different Poisson spike trains for each seed then use:
        # b2.seed(config.randSeed)
        b2.seed(0)
 
        self.network.run(config.duration)
        b2.device.build(
            directory='brian2_standalone_output',
            compile=True,
            run=True,
            clean=True
        )
        print("Simulation complete.\n")

    # =========================================================================
    # ANALYSIS HELPERS
    # =========================================================================

    def _get_spike_arrays(self):
        """
        Return flat (t_ms, gid) arrays for every cell type, ordered to match
        the NEURON gid numbering convention.
        """
        Ns  = config.Npyso
        Nd  = config.Npydr
        Ni  = config.Ninh
        Ntc = config.Ntc
        Nre = config.Nre

        t_all, id_all = [], []

        def _add(mon, offset):
            t_all.extend((mon.t / b2.ms).tolist())
            id_all.extend((np.array(mon.i) + offset).tolist())

        _add(self.PYso_spikemon, 0)
        _add(self.PYdr_spikemon, Ns)
        _add(self.IN_spikemon,   Ns + Nd)
        if self.params['include_thalamus']:
            _add(self.TC_spikemon, Ns + Nd + Ni)
            _add(self.RE_spikemon, Ns + Nd + Ni + Ntc)

        order   = np.argsort(t_all)
        t_arr   = np.array(t_all)[order]
        id_arr  = np.array(id_all)[order]
        return t_arr, id_arr

    # =========================================================================
    # PLOTTING
    # =========================================================================

    def plotRaster(self, output_dir='results'):
        ### Function to create a raster plot after gathering the spikes ###
        print('Plotting raster ...')
        
        # 1. Fetch the flat arrays exactly as NEURON formats them
        tVecAll_arr, idVecAll_arr = self._get_spike_arrays()
        tVecAll = tVecAll_arr.tolist()
        idVecAll = idVecAll_arr.tolist()

        # 2. Local variables for NEURON parity
        Npyso = config.Npyso
        Npydr = config.Npydr
        Ninh  = config.Ninh
        Ntc   = config.Ntc
        Nre   = config.Nre
        N     = Npyso + Npydr + Ninh + Ntc + Nre

        # Indices for each cell type based on their ID ranges
        pyso_indices = [i for i, val in enumerate(idVecAll) if val < Npyso]
        pydr_indices = [i for i, val in enumerate(idVecAll) if Npyso <= val < Npyso + Npydr]
        inh_indices = [i for i, val in enumerate(idVecAll) if Npyso + Npydr <= val < Npyso + Npydr + Ninh]
        if config.current_state_params['include_thalamus']:
            tc_indices = [i for i, val in enumerate(idVecAll) if Npyso + Npydr + Ninh <= val < Npyso + Npydr + Ninh + Ntc]
            re_indices = [i for i, val in enumerate(idVecAll) if Npyso + Npydr + Ninh + Ntc <= val < N]

        # Spike times and IDs for each cell type
        pyso_tVec = [tVecAll[i] for i in pyso_indices]
        pyso_id = [idVecAll[i] for i in pyso_indices]
        pydr_tVec = [tVecAll[i] for i in pydr_indices]
        pydr_id = [idVecAll[i] for i in pydr_indices]
        inh_tVec = [tVecAll[i] for i in inh_indices]
        inh_id = [idVecAll[i] for i in inh_indices]
        if config.current_state_params['include_thalamus']:
            tc_tVec = [tVecAll[i] for i in tc_indices]
            tc_id = [idVecAll[i] for i in tc_indices]
            re_tVec = [tVecAll[i] for i in re_indices]
            re_id = [idVecAll[i] for i in re_indices]

        if config.current_state_params['include_thalamus']:
            # Create subplots: 5 rows (one per cell type), 1 column, share x-axis
            fig, axes = plt.subplots(5, 1, figsize=(9.6, 10.8), sharex=True, constrained_layout=True)
        else:
            # Create subplots: 3 rows (PYso, PYdr, IN), 1 column, share x-axis
            fig, axes = plt.subplots(3, 1, figsize=(9.6, 6.48), sharex=True, constrained_layout=True)

        # List of cell types, their data, colors, and labels
        if config.current_state_params['include_thalamus']:
            cell_types = [
                (pyso_tVec, pyso_id, 'black', 'PYso'),
                (pydr_tVec, pydr_id, 'black', 'PYdr'),
                (inh_tVec, inh_id, 'black', 'IN'),
                (tc_tVec, tc_id, 'black', 'TC'),
                (re_tVec, re_id, 'black', 'RE')
            ]
        else:
            cell_types = [
                (pyso_tVec, pyso_id, 'black', 'PYso'),
                (pydr_tVec, pydr_id, 'black', 'PYdr'),
                (inh_tVec, inh_id, 'black', 'IN')
            ]

        # ID ranges for each cell type (for y-axis limits)
        if config.current_state_params['include_thalamus']:
            id_ranges = [
                (0, Npyso),  # PYso: 0 to Npyso-1
                (Npyso, Npyso + Npydr),  # PYdr: Npyso to Npyso+Npydr-1
                (Npyso + Npydr, Npyso + Npydr + Ninh),  # IN
                (Npyso + Npydr + Ninh, Npyso + Npydr + Ninh + Ntc),  # TC
                (Npyso + Npydr + Ninh + Ntc, N)  # RE
            ]
        else:
            id_ranges = [
                (0, Npyso),  # PYso: 0 to Npyso-1
                (Npyso, Npyso + Npydr),  # PYdr: Npyso to Npyso+Npydr-1
                (Npyso + Npydr, Npyso + Npydr + Ninh)  # IN
            ]

        # Plot each cell type in its subplot
        for ax, (tVec, ids, color, label), (id_min, id_max) in zip(axes, cell_types, id_ranges):
            # Group spike times by cell ID for eventplot (bars)
            spike_dict = {}
            for t, cell_id in zip(tVec, ids):
                t_sec = t / 1000.0
                if cell_id not in spike_dict:
                    spike_dict[cell_id] = []
                spike_dict[cell_id].append(t_sec)

            # Plot bars using eventplot, but only if there are spikes
            spike_times = [spike_dict[cell_id] for cell_id in sorted(spike_dict.keys())]
            cell_ids = sorted(spike_dict.keys())
            
            if spike_times:  # Check if there are any spikes to plot
                # Plot bars
                ax.eventplot(spike_times, lineoffsets=cell_ids, colors=color, linelengths=0.9, linewidths=0.3)
            else:
                # If no spikes, add a text annotation
                ax.text(0.5, 0.5, 'No spikes', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, color='gray')

            # Set y-axis limits and invert (0 at top, max at bottom)
            ax.set_ylim(id_min, id_max)
            ax.invert_yaxis()
            ax.set_ylabel(f'{label}')

        # Set x-axis label on the bottom subplot
        axes[-1].set_xlabel('Time (s)')

        # Set x-axis limits for all subplots
        if tVecAll:
            for ax in axes:
                ax.set_xlim(0, 1.01 * max(tVecAll) / 1000.0)
        
        # --- Save and close ---
        
        os.makedirs(output_dir, exist_ok=True)
        filepath = os.path.join(output_dir, f'raster_plot_{config.CONDITION}.png')
        plt.savefig(filepath, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        plt.close()
        
        print(f'Raster saved: {filepath}')

    def plotVoltageTrace(self, cell_type, cell_idx, output_dir='results'):
        monitors = {
            'PYso': getattr(self, 'PYso_monitor', None),
            'PYdr': getattr(self, 'PYdr_monitor', None),
            'IN':   getattr(self, 'IN_monitor', None),
            'TC':   getattr(self, 'TC_monitor', None),
            'RE':   getattr(self, 'RE_monitor', None)
        }
        
        monitor = monitors.get(cell_type)
        if monitor is None:
            print(f"Warning: No voltage monitor found for {cell_type}.")
            return
            
        t = monitor.t / b2.ms
        v = monitor.v[cell_idx] / b2.mV
        
        plt.figure()
        somaPlot = plt.plot(t, v, color='black')
        plt.legend(somaPlot, ['soma'])
        plt.xlabel('time (ms)')
        plt.ylabel('mV')
        plt.title(f'Cell {cell_type}[{cell_idx}] voltage trace')
        
        os.makedirs(output_dir, exist_ok=True)
        filepath = os.path.join(output_dir, f'traces_cell_{cell_type}_{cell_idx}.png')
        plt.savefig(filepath, bbox_inches='tight')
        plt.close()
        
        print(f'Voltage trace saved: {filepath}')

    # =========================================================================
    # SPIKE RATE REPORTING
    # =========================================================================

    def printSpikeRates(self):
        """Print mean population firing rates for each cell type."""
        sim_s = config.duration / b2.second

        populations = [
            ('PYso', self.PYso_spikemon, config.Npyso),
            ('PYdr', self.PYdr_spikemon, config.Npydr),
            ('IN',   self.IN_spikemon,   config.Ninh),
        ]
        if self.params['include_thalamus']:
            populations += [
                ('TC', self.TC_spikemon, config.Ntc),
                ('RE', self.RE_spikemon, config.Nre),
            ]

        print('\n--- Population Firing Rates ---')
        for name, mon, n_cells in populations:
            n_spikes = len(mon.i)
            rate     = n_spikes / (n_cells * sim_s) if n_cells > 0 else 0.0
            print(f'  {name:5s}: {rate:7.3f} Hz  ({n_spikes} spikes, {n_cells} cells)')
        print()

    # =========================================================================
    # DATA SAVING
    # =========================================================================

    def saveData(self, output_base_dir=os.path.join('results', 'data')):
        print(f'Saving data for condition: {config.CONDITION} ...')

        t_arr, id_arr = self._get_spike_arrays()

        data = {
            'network_params': {
                'Npyso'        : config.Npyso,
                'Npydr'        : config.Npydr,
                'Ninh'         : config.Ninh,
                'Ntc'          : config.Ntc,
                'Nre'          : config.Nre,
                'condition'    : config.CONDITION,
                'stim_type'    : config.stim_type,
                'stim_start'   : config.stim_start,
                'stim_dur'     : config.stim_dur,
                'seed'         : config.randSeed,
                'gAMPA_PY_PY'  : self.params['gAMPA_PYso_PYdr'],
                'gAMPA_TC_PY'  : self.params['gAMPA_TC_PYdr'],
                'recorded_voltages': config.record_voltages,
            },
            'spikes': {
                'tVec' : t_arr.tolist(),
                'idVec': id_arr.tolist(),
            },
        }

        if config.record_voltages:
            all_voltages = {}
            Ns = config.Npyso
            Nd = config.Npydr

            for idx in range(Ns):
                all_voltages[idx] = self.PYso_monitor.v[idx] / b2.mV

            eeg_signal = np.zeros(len(self.PYdr_monitor.t))
            for idx in range(Nd):
                v_trace = self.PYdr_monitor.v[idx] / b2.mV
                all_voltages[Ns + idx] = v_trace
                eeg_signal += v_trace

            data['voltages'] = all_voltages
            data['eeg'] = {
                'time'  : self.PYdr_monitor.t / b2.ms,
                'signal': eeg_signal,
            }
            print(' -> Included full voltage arrays and computed EEG proxy.')
        else:
            print(' -> Included minimal raster data only.')

        if config.CONDITION in ('propofol_mostly_c', 'propofol_mostly_i') \
                and config.stim_type == 'none':
            condition_dir = output_base_dir
        else:
            condition_dir = os.path.join(output_base_dir, config.CONDITION)
        os.makedirs(condition_dir, exist_ok=True)

        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        if config.args.table_name in ('table1', 'table2'):
            gTC = f"{self.params['gAMPA_TC_PYdr']*1e3:.4f}"
            gPY = f"{self.params['gAMPA_PYso_PYdr']*1e3:.4f}"
            filename = (f"sim_{config.args.table_name}"
                        f"_gTC_{gTC}_gPY_{gPY}_seed_{config.randSeed}.pkl")
        elif config.CONDITION in ('propofol_mostly_c', 'propofol_mostly_i'):
            filename = (f"sim_{config.CONDITION}_{config.stim_type}"
                        f"_{int(config.stim_start)}_{timestamp}.pkl")
        else:
            filename = f'sim_data_{config.CONDITION}_{timestamp}.pkl'

        filepath = os.path.join(condition_dir, filename)
        with open(filepath, 'wb') as fh:
            pickle.dump(data, fh)
        print(f'Data saved → {filepath}')
        return filepath
