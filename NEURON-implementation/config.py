'''
Define baseline parameters for model with condition switching.
'''

from neuron import h
# from neuron import gui # enable when not running on cluster
import numpy as np
import argparse

# ============================================================================
# COMMAND LINE CONFIGURATION
# ============================================================================
parser = argparse.ArgumentParser(description='Propofol Network Simulation')
parser.add_argument('--condition', type=str, default='propofol', choices=['wake', 'cortex_only_swo', 'propofol', 'disconnected_thalamus',
                                                                                    'propofol_mostly_c', 'propofol_mostly_i', 'tables_propofol'])
parser.add_argument('--stim_type', type=str, default='none', choices=['none', 'sync', 'desync'])
parser.add_argument('--stim_start', type=float, default=7500.0) # ms
parser.add_argument('--record_voltages', action='store_true', help='Flag to record continuous voltages')

# --- NOTE ON --record_voltages FLAG ---
# This will save the full voltage traces for all cells, which can be very large. 
# Use with caution, especially for long simulations or large networks. The default is False to save disk space and speed up simulations. Also, 
# when running with --record_voltages, please do not use --coreneuron or --gpu flags, as voltage recording is not currently supported in those modes.
# Insted, run with the standard NEURON backend to ensure voltage recording works correctly.
# Compile all the mod files with standard NEURON before running with --record_voltages, and do not use CoreNEURON or GPU mode for those runs.

parser.add_argument('--coreneuron', action='store_true', help='Flag to run with CoreNEURON')
parser.add_argument('--gpu', action='store_true', help='Flag to run CoreNEURON on GPU (requires GPU build)')
# parser.add_argument('--duration', type=float, default=30000.0) # ms

# Arguments for Tables 1 and 2
parser.add_argument('--table_name', type=str, default='none', help='Identifier for table batch runs (e.g., table1, table2)')
parser.add_argument('--gAMPA_PY_PY', type=float, default=None, help='Overrides gAMPA_PYso_PYdr')
parser.add_argument('--gAMPA_TC_PY', type=float, default=None, help='Overrides gAMPA_TC_PYdr')
parser.add_argument('--gext', type=float, default=None, help='Overrides background Poisson noise gext_PYdr, gext_TC, and gext_RE')
parser.add_argument('--seed', type=int, default=101, help='Random seed for independent trials')

# Parse known args so NEURON's internal MPI args don't cause crashes
args, unknown = parser.parse_known_args()

CONDITION = args.condition
stim_type = args.stim_type
stim_start = args.stim_start
stim_dur = 100.0 # ms
record_voltages = args.record_voltages

# Override the global seed with the command line argument
randSeed = args.seed # randomizer seed for random number generation
# ============================================================================

# --- ParallelContext object --- 
pc = h.ParallelContext()
pc.set_maxstep(10)  # ms, maximum step size for parallel simulation
idhost = int(pc.id())
nhost = int(pc.nhost())

# duration of simulation
duration = 30000 # (ms)
t_seg = 10.0     # (ms) simulation time between each data dump to node 0
time_step = 0.01 # (ms) time step for simulation

# # set randomizer seed
# randSeed = 101 # global seed for random number generation

# set numbers of each cell type
Npydr = 100
Npyso = 100
Ninh = 20
Nre = 20
Ntc = 20

# threshold for detecting spikes
threshold = -25

# --- Poisson simulation parameters ---
dt = h.dt  # ms, NEURON time step

# Common parameters across all conditions
# For PYdr cells
baseline_PYdr = 40.0  # Hz
DC_PYdr = 0           # Hz 
AC_PYdr = 0           # Hz 
f_PYdr = 0            # Hz 
phi_PYdr = 0          # radians
onset_PYdr = [0]      # ms
offset_PYdr = [np.inf]# ms
tau_PYdr = 2          # ms
Eext_PYdr = 0         # mV (reversal potential)

# For TC cells
baseline_TC = 40.0    # Hz
DC_TC = 0             # Hz
AC_TC = 0             # Hz
f_TC = 0              # Hz
phi_TC = 0            # radians
onset_TC = [0]        # ms
offset_TC = [np.inf]  # ms
tau_TC = 2            # ms
Eext_TC = 0           # mV (reversal potential)

# For RE cells
baseline_RE = 40.0    # Hz
DC_RE = 0             # Hz
AC_RE = 0             # Hz
f_RE = 0              # Hz
phi_RE = 0            # radians
onset_RE = [0]        # ms
offset_RE = [np.inf]  # ms
tau_RE = 2            # ms
Eext_RE = 0           # mV (reversal potential)

# === STP parameters ===
deprFactor = 0.9      # depression factor (0 < deprFactor <= 1)
tauRes = 400          # recovery time constant (ms)
release_time = 2*dt   # delay for depression effect (ms)

# ============================================================================
# CONDITION-SPECIFIC PARAMETERS
# ============================================================================

class ConditionParams:
    """Store parameters for each condition"""
    
    def __init__(self):
        # Define all the conditions
        self.conditions = {
            'wake': {
                'gext_PYdr': 0.004e-3,
                'gext_TC': 0.004e-3,
                'gext_RE': 0.004e-3,
                'gH': 0.4e-3,
                'gKLeak': 0.001e-3,
                'gKNa': 0.1e-3,
                'gAMPA_PYso_PYdr': 0.004e-3,
                'gNMDA_PYso_PYdr': 0.00257e-3 * 0.4,
                'gAMPA_PYso_IN': 1.0e-3 * 0.4,
                'gNMDA_PYso_IN': 0.0025e-3 * 0.4,
                'gGABAA_IN_PYso': 0.1e-3 * 0.4,
                'gAMPA_TC_PYdr': 0.004e-3,
                'gAMPA_TC_IN': 0.1e-3 * 0.4,
                'gAMPA_PYso_TC': 0.4e-3 * 0.4,
                'gAMPA_PYso_RE': 0.2e-3 * 0.4,
                'propoCondMult': 1,
                'propoTauMult': 1,
                'include_thalamus': True,
                'include_thalamocortical_connections': True
            },
            'cortex_only_swo': {
                'gext_PYdr': 0.005e-3,
                'gext_TC': 0.005e-3,
                'gext_RE': 0.005e-3,
                'gH': 0.005e-3,
                'gKLeak': 0.0172e-3,
                'gKNa': 1.33e-3,
                'gAMPA_PYso_PYdr': 0.005e-3,
                'gNMDA_PYso_PYdr': 0.00257e-3,
                'gAMPA_PYso_IN': 1.0e-3,
                'gNMDA_PYso_IN': 0.0025e-3,
                'gGABAA_IN_PYso': 0.1e-3,
                'gAMPA_TC_PYdr': 0.005e-3,
                'gAMPA_TC_IN': 0.1e-3,
                'gAMPA_PYso_TC': 0.4e-3,
                'gAMPA_PYso_RE': 0.2e-3,
                'propoCondMult': 1,
                'propoTauMult': 1,
                'include_thalamus': False,  # <--- This is the key difference
                'include_thalamocortical_connections': False
            },
            'propofol': {
                'gext_PYdr': 0.005e-3,
                'gext_TC': 0.005e-3,
                'gext_RE': 0.005e-3,
                'gH': 0.005e-3,
                'gKLeak': 0.0172e-3,
                'gKNa': 1.33e-3,
                'gAMPA_PYso_PYdr': 0.005e-3,
                'gNMDA_PYso_PYdr': 0.00257e-3,
                'gAMPA_PYso_IN': 1.0e-3,
                'gNMDA_PYso_IN': 0.0025e-3,
                'gGABAA_IN_PYso': 0.1e-3,
                'gAMPA_TC_PYdr': 0.005e-3,
                'gAMPA_TC_IN': 0.1e-3,
                'gAMPA_PYso_TC': 0.4e-3,
                'gAMPA_PYso_RE': 0.2e-3,
                'propoCondMult': 3,
                'propoTauMult': 3,
                'include_thalamus': True,
                'include_thalamocortical_connections': True
            },
            'disconnected_thalamus': {
                'gext_PYdr': 0.005e-3,
                'gext_TC': 0.005e-3,
                'gext_RE': 0.005e-3,
                'gH': 0.005e-3,
                'gKLeak': 0.0172e-3,
                'gKNa': 1.33e-3,
                'gAMPA_PYso_PYdr': 0.008e-3,
                'gNMDA_PYso_PYdr': 0.00257e-3,
                'gAMPA_PYso_IN': 1.0e-3,
                'gNMDA_PYso_IN': 0.0025e-3,
                'gGABAA_IN_PYso': 0.1e-3,
                'gAMPA_TC_PYdr': 0,  # Disconnected thalamus
                'gAMPA_TC_IN': 0,    # Disconnected thalamus
                'gAMPA_PYso_TC': 0,  # Disconnected thalamus
                'gAMPA_PYso_RE': 0,  # Disconnected thalamus
                'propoCondMult': 3,
                'propoTauMult': 3,
                'include_thalamus': True,   # Thalamus is present but disconnected
                'include_thalamocortical_connections': False
            },
            'propofol_mostly_c': {
                'gext_PYdr': 0.005e-3, 'gext_TC': 0.005e-3, 'gext_RE': 0.005e-3,
                'gH': 0.005e-3, 'gKLeak': 0.0172e-3, 'gKNa': 1.33e-3,
                'gAMPA_PYso_PYdr': 0.002e-3,  # Lowered to favor C-state
                'gNMDA_PYso_PYdr': 0.00257e-3,
                'gAMPA_PYso_IN': 1.0e-3, 'gNMDA_PYso_IN': 0.0025e-3, 'gGABAA_IN_PYso': 0.1e-3,
                'gAMPA_TC_PYdr': 0.002e-3,    # Lowered to favor C-state
                'gAMPA_TC_IN': 0.1e-3, 'gAMPA_PYso_TC': 0.4e-3, 'gAMPA_PYso_RE': 0.2e-3,
                'propoCondMult': 3, 'propoTauMult': 3,
                'include_thalamus': True, 'include_thalamocortical_connections': True
            },
            'propofol_mostly_i': {
                'gext_PYdr': 0.005e-3, 'gext_TC': 0.005e-3, 'gext_RE': 0.005e-3,
                'gH': 0.005e-3, 'gKLeak': 0.0172e-3, 'gKNa': 1.33e-3,
                'gAMPA_PYso_PYdr': 0.008e-3,  # Raised to force I-state
                'gNMDA_PYso_PYdr': 0.00257e-3,
                'gAMPA_PYso_IN': 1.0e-3, 'gNMDA_PYso_IN': 0.0025e-3, 'gGABAA_IN_PYso': 0.1e-3,
                'gAMPA_TC_PYdr': 0.008e-3,    # Raised to force I-state
                'gAMPA_TC_IN': 0.1e-3, 'gAMPA_PYso_TC': 0.4e-3, 'gAMPA_PYso_RE': 0.2e-3,
                'propoCondMult': 3, 'propoTauMult': 3,
                'include_thalamus': True, 'include_thalamocortical_connections': True
            }
        }
    
    def get_params(self, condition):
        """Get parameters for specified condition"""
        if condition not in self.conditions:
            raise ValueError(f"Unknown condition: {condition}. Available: {list(self.conditions.keys())}")
        return self.conditions[condition]

# Create the condition parameters object
_condition_params = ConditionParams()

# ============================================================================
# PARAMETER FETCHING & OVERRIDES
# ============================================================================

# Create an alias for the tables batch runs to route them to a separate folder
if CONDITION == 'tables_propofol':
    current_state_params = _condition_params.get_params('propofol').copy()
else:
    current_state_params = _condition_params.get_params(CONDITION)

# --- Override specific conductances if provided via command line ---
if args.gAMPA_PY_PY is not None:
    current_state_params['gAMPA_PYso_PYdr'] = args.gAMPA_PY_PY
if args.gAMPA_TC_PY is not None:
    current_state_params['gAMPA_TC_PYdr'] = args.gAMPA_TC_PY
# ------------------------------------------------------------------------

# --- Override the Poisson noise conductances if --gext is provided ---
if args.gext is not None:
    current_state_params['gext_PYdr'] = args.gext
    current_state_params['gext_TC'] = args.gext
    current_state_params['gext_RE'] = args.gext
# ------------------------------------------------------------------------

# Print current condition on rank 0
if idhost == 0:
    print("=" * 70)
    print(f"SIMULATION CONDITION: {CONDITION.upper()}")
    print("=" * 70)
    print(f"Include thalamus: {current_state_params['include_thalamus']}")
    print(f"Include thalamocortical connections: {current_state_params['include_thalamocortical_connections']}")
    print(f"gext_PYdr: {current_state_params['gext_PYdr']:.6e}")
    print(f"gH: {current_state_params['gH']:.6e}")
    print(f"gKLeak: {current_state_params['gKLeak']:.6e}")
    print(f"propoCondMult: {current_state_params['propoCondMult']}")
    print(f"propoTauMult: {current_state_params['propoTauMult']}")
    print("=" * 70)
