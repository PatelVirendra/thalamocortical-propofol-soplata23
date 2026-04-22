"""
config.py — Brian2 Thalamocortical Propofol Model
==================================================
Defines all global simulation parameters, condition-specific conductance
tables, CLI argument parsing, and Brian2 device/clock setup.

Mirrors the structure of the NEURON config.py so that parameter names,
condition keys, and CLI flags are identical between the two implementations.

Reference: Soplata et al. (2023)
"""

import argparse
import numpy as np
import brian2 as b2

# =============================================================================
# COMMAND-LINE ARGUMENT PARSING
# =============================================================================

parser = argparse.ArgumentParser(description='Brian2 Propofol Thalamocortical Network Simulation')

parser.add_argument(
    '--condition', type=str, default='propofol',
    choices=['wake', 'cortex_only_swo', 'propofol', 'disconnected_thalamus',
             'propofol_mostly_c', 'propofol_mostly_i', 'tables_propofol'],
    help='Simulation condition (determines conductance parameters)'
)
parser.add_argument(
    '--stim_type', type=str, default='none', choices=['none', 'sync', 'desync'],
    help='Type of external stimulation applied during simulation'
)
parser.add_argument(
    '--stim_start', type=float, default=7500.0,
    help='Stimulation onset time (ms)'
)
parser.add_argument(
    '--record_voltages', action='store_true',
    help='If set, record continuous voltage traces for all cells (large output)'
)
parser.add_argument(
    '--seed', type=int, default=101,
    help='Random seed for reproducibility across independent trials'
)
parser.add_argument(
    '--openmp_threads', type=int, default=2,
    help='Number of OpenMP threads for C++ standalone compilation'
)

# --- Arguments for Tables 1 & 2 batch runs ---
parser.add_argument(
    '--table_name', type=str, default='none',
    help='Identifier for table batch runs (e.g., table1, table2)'
)
parser.add_argument(
    '--gAMPA_PY_PY', type=float, default=None,
    help='Override gAMPA_PYso_PYdr (S/cm²)'
)
parser.add_argument(
    '--gAMPA_TC_PY', type=float, default=None,
    help='Override gAMPA_TC_PYdr (S/cm²)'
)
parser.add_argument(
    '--gext', type=float, default=None,
    help='Override background Poisson noise conductance for PYdr, TC, and RE'
)

args, _unknown = parser.parse_known_args()

# =============================================================================
# TOP-LEVEL SIMULATION SETTINGS
# =============================================================================

CONDITION      = args.condition
stim_type      = args.stim_type
stim_start     = args.stim_start   # ms
stim_dur       = 100.0             # ms  (fixed, not exposed on CLI)
record_voltages = args.record_voltages
randSeed       = args.seed

# --- Simulation timing ---
duration        = 30000 * b2.ms    # total run time
time_step       = 0.01  * b2.ms   # integration time step (0.01 ms = 10 kHz)
downsample_factor = 10            # record every N-th step → 1 kHz output
record_interval   = downsample_factor * time_step

# --- Population sizes ---
Npydr = 100
Npyso = 100
Ninh  = 20
Nre   = 20
Ntc   = 20

# --- Spike detection threshold ---
threshold_mV = -25   # mV

# --- Synaptic transmission delays ---
transmission_delay = 0.01 * b2.ms   # axonal transmission delay
release_time       = 2 * time_step  # STP depression window width

# =============================================================================
# BACKGROUND POISSON NOISE PARAMETERS  (shared across all conditions)
# =============================================================================

# PYdr background drive
baseline_PYdr = 40.0 * b2.Hz
tau_PYdr      = 2    * b2.ms
Eext_PYdr     = 0    * b2.mV

# TC background drive
baseline_TC   = 40.0 * b2.Hz
tau_TC        = 2    * b2.ms
Eext_TC       = 0    * b2.mV

# RE background drive
baseline_RE   = 40.0 * b2.Hz
tau_RE        = 2    * b2.ms
Eext_RE       = 0    * b2.mV

# =============================================================================
# CONDITION-SPECIFIC CONDUCTANCE PARAMETERS
# =============================================================================

class ConditionParams:
    """
    Stores per-condition conductance tables.

    Keys match NEURON config.py exactly so downstream analysis scripts
    (HMM classifiers, table-generation scripts, etc.) are condition-agnostic.
    All raw conductances are in S/cm².
    """

    def __init__(self):
        self.conditions = {

            # ----------------------------------------------------------------
            'wake': {
                'gext_PYdr'          : 0.004e-3,
                'gext_TC'            : 0.004e-3,
                'gext_RE'            : 0.004e-3,
                'gH'                 : 0.4e-3,
                'gKLeak'             : 0.001e-3,
                'gKNa'               : 0.1e-3,
                'gAMPA_PYso_PYdr'    : 0.004e-3,
                'gNMDA_PYso_PYdr'    : 0.00257e-3 * 0.4,
                'gAMPA_PYso_IN'      : 1.0e-3 * 0.4,
                'gNMDA_PYso_IN'      : 0.0025e-3 * 0.4,
                'gGABAA_IN_PYso'     : 0.1e-3 * 0.4,
                'gAMPA_TC_PYdr'      : 0.004e-3,
                'gAMPA_TC_IN'        : 0.1e-3 * 0.4,
                'gAMPA_PYso_TC'      : 0.4e-3 * 0.4,
                'gAMPA_PYso_RE'      : 0.2e-3 * 0.4,
                'propoCondMult'      : 1,
                'propoTauMult'       : 1,
                'include_thalamus'               : True,
                'include_thalamocortical_connections' : True,
            },

            # ----------------------------------------------------------------
            'cortex_only_swo': {
                'gext_PYdr'          : 0.005e-3,
                'gext_TC'            : 0.005e-3,
                'gext_RE'            : 0.005e-3,
                'gH'                 : 0.005e-3,
                'gKLeak'             : 0.0172e-3,
                'gKNa'               : 1.33e-3,
                'gAMPA_PYso_PYdr'    : 0.005e-3,
                'gNMDA_PYso_PYdr'    : 0.00257e-3,
                'gAMPA_PYso_IN'      : 1.0e-3,
                'gNMDA_PYso_IN'      : 0.0025e-3,
                'gGABAA_IN_PYso'     : 0.1e-3,
                'gAMPA_TC_PYdr'      : 0.005e-3,
                'gAMPA_TC_IN'        : 0.1e-3,
                'gAMPA_PYso_TC'      : 0.4e-3,
                'gAMPA_PYso_RE'      : 0.2e-3,
                'propoCondMult'      : 1,
                'propoTauMult'       : 1,
                'include_thalamus'               : False,   # cortex-only
                'include_thalamocortical_connections' : False,
            },

            # ----------------------------------------------------------------
            'propofol': {
                'gext_PYdr'          : 0.005e-3,
                'gext_TC'            : 0.005e-3,
                'gext_RE'            : 0.005e-3,
                'gH'                 : 0.005e-3,
                'gKLeak'             : 0.0172e-3,
                'gKNa'               : 1.33e-3,
                'gAMPA_PYso_PYdr'    : 0.005e-3,
                'gNMDA_PYso_PYdr'    : 0.00257e-3,
                'gAMPA_PYso_IN'      : 1.0e-3,
                'gNMDA_PYso_IN'      : 0.0025e-3,
                'gGABAA_IN_PYso'     : 0.1e-3,
                'gAMPA_TC_PYdr'      : 0.005e-3,
                'gAMPA_TC_IN'        : 0.1e-3,
                'gAMPA_PYso_TC'      : 0.4e-3,
                'gAMPA_PYso_RE'      : 0.2e-3,
                'propoCondMult'      : 3,
                'propoTauMult'       : 3,
                'include_thalamus'               : True,
                'include_thalamocortical_connections' : True,
            },

            # ----------------------------------------------------------------
            'disconnected_thalamus': {
                'gext_PYdr'          : 0.005e-3,
                'gext_TC'            : 0.005e-3,
                'gext_RE'            : 0.005e-3,
                'gH'                 : 0.005e-3,
                'gKLeak'             : 0.0172e-3,
                'gKNa'               : 1.33e-3,
                'gAMPA_PYso_PYdr'    : 0.008e-3,
                'gNMDA_PYso_PYdr'    : 0.00257e-3,
                'gAMPA_PYso_IN'      : 1.0e-3,
                'gNMDA_PYso_IN'      : 0.0025e-3,
                'gGABAA_IN_PYso'     : 0.1e-3,
                'gAMPA_TC_PYdr'      : 0.0,   # thalamus present but silent
                'gAMPA_TC_IN'        : 0.0,
                'gAMPA_PYso_TC'      : 0.0,
                'gAMPA_PYso_RE'      : 0.0,
                'propoCondMult'      : 3,
                'propoTauMult'       : 3,
                'include_thalamus'               : True,   # cells exist
                'include_thalamocortical_connections' : False,  # but disconnected
            },

            # ----------------------------------------------------------------
            'propofol_mostly_c': {
                'gext_PYdr'          : 0.005e-3,
                'gext_TC'            : 0.005e-3,
                'gext_RE'            : 0.005e-3,
                'gH'                 : 0.005e-3,
                'gKLeak'             : 0.0172e-3,
                'gKNa'               : 1.33e-3,
                'gAMPA_PYso_PYdr'    : 0.002e-3,   # lowered → favours C-state
                'gNMDA_PYso_PYdr'    : 0.00257e-3,
                'gAMPA_PYso_IN'      : 1.0e-3,
                'gNMDA_PYso_IN'      : 0.0025e-3,
                'gGABAA_IN_PYso'     : 0.1e-3,
                'gAMPA_TC_PYdr'      : 0.002e-3,   # lowered → favours C-state
                'gAMPA_TC_IN'        : 0.1e-3,
                'gAMPA_PYso_TC'      : 0.4e-3,
                'gAMPA_PYso_RE'      : 0.2e-3,
                'propoCondMult'      : 3,
                'propoTauMult'       : 3,
                'include_thalamus'               : True,
                'include_thalamocortical_connections' : True,
            },

            # ----------------------------------------------------------------
            'propofol_mostly_i': {
                'gext_PYdr'          : 0.005e-3,
                'gext_TC'            : 0.005e-3,
                'gext_RE'            : 0.005e-3,
                'gH'                 : 0.005e-3,
                'gKLeak'             : 0.0172e-3,
                'gKNa'               : 1.33e-3,
                'gAMPA_PYso_PYdr'    : 0.008e-3,   # raised → forces I-state
                'gNMDA_PYso_PYdr'    : 0.00257e-3,
                'gAMPA_PYso_IN'      : 1.0e-3,
                'gNMDA_PYso_IN'      : 0.0025e-3,
                'gGABAA_IN_PYso'     : 0.1e-3,
                'gAMPA_TC_PYdr'      : 0.008e-3,   # raised → forces I-state
                'gAMPA_TC_IN'        : 0.1e-3,
                'gAMPA_PYso_TC'      : 0.4e-3,
                'gAMPA_PYso_RE'      : 0.2e-3,
                'propoCondMult'      : 3,
                'propoTauMult'       : 3,
                'include_thalamus'               : True,
                'include_thalamocortical_connections' : True,
            },
        }

    def get_params(self, condition):
        if condition not in self.conditions:
            raise ValueError(
                f"Unknown condition: '{condition}'. "
                f"Available: {list(self.conditions.keys())}"
            )
        return self.conditions[condition]


# =============================================================================
# ACTIVE PARAMETER RESOLUTION  (with CLI overrides)
# =============================================================================

_condition_params = ConditionParams()

# 'tables_propofol' reuses the 'propofol' parameters but saves to a separate folder
if CONDITION == 'tables_propofol':
    current_state_params = _condition_params.get_params('propofol').copy()
else:
    current_state_params = _condition_params.get_params(CONDITION).copy()

# --- Per-parameter CLI overrides (for batch sweeps, Tables 1 & 2) ---
if args.gAMPA_PY_PY is not None:
    current_state_params['gAMPA_PYso_PYdr'] = args.gAMPA_PY_PY
if args.gAMPA_TC_PY is not None:
    current_state_params['gAMPA_TC_PYdr'] = args.gAMPA_TC_PY
if args.gext is not None:
    current_state_params['gext_PYdr'] = args.gext
    current_state_params['gext_TC']   = args.gext
    current_state_params['gext_RE']   = args.gext

# =============================================================================
# BRIAN2 DEVICE & CLOCK SETUP
# =============================================================================

# Reinitialise to allow multiple imports (e.g. during testing)
b2.device.reinit()
b2.device.activate(
    directory='brian2_standalone_output',
    build_on_run=False
)
b2.set_device(
    'cpp_standalone',
    directory='brian2_standalone_output',
    build_on_run=False
)
b2.prefs.devices.cpp_standalone.openmp_threads = args.openmp_threads
b2.defaultclock.dt = time_step

# Convenience unit for micromolar concentrations used in several cell equations
uM = 1e-3 * b2.mM

# =============================================================================
# STARTUP BANNER
# =============================================================================

print("=" * 70)
print(f"SIMULATION CONDITION : {CONDITION.upper()}")
print("=" * 70)
print(f"Include thalamus                : {current_state_params['include_thalamus']}")
print(f"Include TC connections          : {current_state_params['include_thalamocortical_connections']}")
print(f"gext_PYdr (S/cm²)               : {current_state_params['gext_PYdr']:.6e}")
print(f"gH (S/cm²)                      : {current_state_params['gH']:.6e}")
print(f"gKLeak (S/cm²)                  : {current_state_params['gKLeak']:.6e}")
print(f"gKNa (S/cm²)                    : {current_state_params['gKNa']:.6e}")
print(f"gAMPA_PYso_PYdr (S/cm²)         : {current_state_params['gAMPA_PYso_PYdr']:.6e}")
print(f"gAMPA_TC_PYdr (S/cm²)           : {current_state_params['gAMPA_TC_PYdr']:.6e}")
print(f"propoCondMult                   : {current_state_params['propoCondMult']}")
print(f"propoTauMult                    : {current_state_params['propoTauMult']}")
print(f"Random seed                     : {randSeed}")
print(f"OpenMP threads                  : {args.openmp_threads}")
print("=" * 70)
