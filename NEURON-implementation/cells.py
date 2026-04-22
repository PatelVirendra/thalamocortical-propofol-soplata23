import config
from neuron import h
import numpy as np
import matplotlib.pyplot as plt
h('load_file("stdgui.hoc")')

class Cell(object):
    """Base class for all neuron types with only a soma section."""
    def __init__(self, gid):
        self.gid = gid
        self.gid_offset = config.Npydr + config.Npyso + config.Ninh + config.Nre + config.Ntc
        self.synlist = []  # List to store synapse objects
        self.createSections()
        self.defineGeometry()
        self.defineBiophysics()
        self.applyPoisson()
        self.applyStimulus()
        self.createSynapses()
        self.createInDegVars()
        self.nclist = []  # List for NetCon objects (if needed later)

    def createSections(self):
        """Create only the soma section."""
        self.soma = h.Section(name='soma', cell=self)

    def defineGeometry(self):
        """Define default geometry for the soma (overridden by subclasses if needed)."""
        pass

    def defineBiophysics(self):
        """Define biophysical properties (overridden by subclasses)."""
        pass

    def applyPoisson(self):
        """Apply background Poisson noise (empty by default, can be overridden)."""
        pass

    def applyStimulus(self):
        """Apply external current stimulus (empty by default, can be overridden)."""
        pass

    def createSynapses(self):
        """Create synapses (overridden by subclasses)."""
        pass

    def createInDegVars(self):
        """Initialize variables to track incoming connections (overridden by subclasses)."""
        pass

    def associateGid(self):
        """Assign a unique gid using parallel context."""
        config.pc.set_gid2node(self.gid, config.idhost)
        config.pc.set_gid2node(self.gid + self.gid_offset, config.idhost)
        nc = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        nc.threshold = config.threshold
        config.pc.cell(self.gid, nc)
        del nc  # NetCon is discarded after registration

    def createNetcon(self, thresh=0):
        """Create a NetCon to record spikes from the soma."""
        nc = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        nc.threshold = config.threshold
        return nc

    def setRecording(self):
        # print(f"Setting up recording for cell {self.gid}")
        """Set up voltage and time recording vectors for the soma."""
        self.soma_v_vec = h.Vector()
        self.tVec = h.Vector()
        self.soma_v_vec.record(self.soma(0.5)._ref_v, 0.1) # Add 0.1 ms sampling!
        self.tVec.record(h._ref_t, 0.1)                    # Add 0.1 ms sampling!

    def plotTraces(self):
        # Convert NEURON Vectors to NumPy arrays
        t = self.tVec.as_numpy()
        v = self.soma_v_vec.as_numpy()
        
        # Create the plot
        plt.figure()
        somaPlot = plt.plot(t, v, color='black')
        plt.legend(somaPlot, ['soma'])
        plt.xlabel('time (ms)')
        plt.ylabel('mV')
        plt.title(f'Cell {self.gid} voltage trace')
        plt.savefig(f'traces_cell_{self.gid}.png')
        plt.close()

    def createIClamp(self, amp):
        """Add a current clamp or an external electrical stimulus to the soma."""
        self.Iext = h.IClamp(self.soma(0.5))
        self.Iext.dur = 1e9  # Effectively infinite duration
        self.Iext.delay = 0
        self.Iext.amp = amp

    def __str__(self):
        return f'Cell[{self.gid}]'


class PYsoCell(Cell):
    def defineGeometry(self):
        """Set soma geometry (default values since not specified)."""
        self.soma.diam = self.soma.L = 69.1   # Length and diameter in microns
        self.soma.nseg = 1

    def defineBiophysics(self):
        """Insert the PYso mechanism into the soma."""
        self.soma.insert('PYso')  # Mechanism defined in PYso.mod

    def createSynapses(self):
        pass

    def createInDegVars(self):
        """Track incoming connections."""
        pass

class PYdrCell(Cell):
    def defineGeometry(self):
        self.soma.L = self.soma.diam = 105.5 # for 0.00035 cm² area
        self.soma.nseg = 1

    def defineBiophysics(self):
        self.soma.insert('PYdr')

    def applyPoisson(self):
        """Applies CVode-safe Poisson noise using NetStim and ExpSyn."""
        
        # 1. Set parameters for the noise from your config file
        rate = config.baseline_PYdr  # e.g., 40 Hz
        
        # 2. Create the NetStim spike generator
        # We attach it to `self` to prevent it from being garbage collected
        self.ns = h.NetStim()
        self.ns.interval = 1000.0 / rate  # Convert rate (Hz) to mean interval (ms)
        self.ns.number = 1e9                     # Large number for continuous firing
        self.ns.noise = 1                 # A value of 1 makes the intervals exponential (a Poisson process)
        self.ns.start = 0                 # Start generating spikes at t=0

        # 3. Create a simple exponential synapse to receive the noise
        # ExpSyn is a built-in NEURON mechanism that does exactly what your Python code did.
        self.syn = h.ExpSyn(self.soma(0.5))
        self.syn.e = config.Eext_PYdr    # Reversal potential for the synapse
        
        # This is the time constant of the conductance decay.
        # It should be the SAME as the `tau` you used in your python script.
        self.syn.tau = config.tau_PYdr   # e.g., 2 ms

        # 4. Connect the spike generator to the synapse with a NetCon
        self.nc = h.NetCon(self.ns, self.syn)
        self.nc.delay = 0.1  # Use a minimal, non-zero delay for numerical stability

        # The 'weight' of the NetCon corresponds to the peak conductance increase (g_max) 
        # for a single noisy event. This is the equivalent of the 'kick' in your python script.
        # We will calculate a good starting value from your old `gext_PYdr` parameter.

        area_um2 = self.soma.L * np.pi * self.soma.diam
        area_cm2 = area_um2 * 1e-8
        
        # Convert gext (S/cm^2) to an absolute conductance weight in microsiemens (uS)
        # The final scaling factor may need to be adjusted to match the old model's firing rate.
        # The term `rate * tau * 1e-3` approximates the average conductance.
        # Let's start with a value derived from your gext parameter.
        g_Abs_uS = config.current_state_params['gext_PYdr'] * area_cm2 * 1e6 

        # The weight needs to be tuned. A good heuristic is to scale the absolute conductance
        # by a factor related to the average number of simultaneous events.
        # Let's start simple: the kick value.
        kick = 1.0 # From your getPoissonGating call
        
        # A reasonable weight is a fraction of the total absolute conductance. Let's start with a small value.
        # The exact value needs to be tuned to match the behavior of your old model.
        # A good starting point for the NetCon weight (in uS) is often a small number like 0.005 to 0.05.
        weight_uS = g_Abs_uS # (uS) <--- THIS IS THE KEY PARAMETER TO TUNE
        
        self.nc.weight[0] = weight_uS

        # print(f"Attached NetStim-based Poisson noise to PYdr cell {self.gid}")

    def applyStimulus(self):
        """Apply external current stimulus to the PYdr cell."""
        pass

    def createSynapses(self):
        pass

    def createInDegVars(self):
        pass

class INCell(Cell):
    def defineGeometry(self):
        self.soma.diam = self.soma.L = 80   # Length and diameter in microns
        self.soma.nseg = 1

    def defineBiophysics(self):
        self.soma.insert('IN')

    def createSynapses(self):
        pass

    def createInDegVars(self):
        pass


class TCCell(Cell):
    def defineGeometry(self):
        self.soma.L = self.soma.diam = np.sqrt(2.9e4/np.pi) # one compartment of 29,000 um2 
        self.soma.nseg = 1

    def defineBiophysics(self):
        self.soma.insert('TC')
        self.soma(0.5).gH_TC = config.current_state_params['gH'] # S/cm2
        self.soma(0.5).gKLeak_TC = config.current_state_params['gKLeak'] # S/cm2

    def applyPoisson(self):
        """Applies CVode-safe Poisson noise using NetStim and ExpSyn."""
        
        # 1. Set parameters for the noise from your config file
        rate = config.baseline_TC  # e.g., 40 Hz
        
        # 2. Create the NetStim spike generator
        # We attach it to `self` to prevent it from being garbage collected
        self.ns = h.NetStim()
        self.ns.interval = 1000.0 / rate  # Convert rate (Hz) to mean interval (ms)
        self.ns.number = 1e9                     # Large number for continuous firing
        self.ns.noise = 1                 # A value of 1 makes the intervals exponential (a Poisson process)
        self.ns.start = 0                 # Start generating spikes at t=0

        # 3. Create a simple exponential synapse to receive the noise
        # ExpSyn is a built-in NEURON mechanism that does exactly what your Python code did.
        self.syn = h.ExpSyn(self.soma(0.5))
        self.syn.e = config.Eext_TC   # Reversal potential for the synapse
        
        # This is the time constant of the conductance decay.
        # It should be the SAME as the `tau` you used in your python script.
        self.syn.tau = config.tau_TC   # e.g., 2 ms

        # 4. Connect the spike generator to the synapse with a NetCon
        self.nc = h.NetCon(self.ns, self.syn)
        self.nc.delay = 0.1  # Use a minimal, non-zero delay for numerical stability

        # The 'weight' of the NetCon corresponds to the peak conductance increase (g_max) 
        # for a single noisy event. This is the equivalent of the 'kick' in your python script.
        # We will calculate a good starting value from your old `gext_PYdr` parameter.

        area_um2 = self.soma.L * np.pi * self.soma.diam
        area_cm2 = area_um2 * 1e-8
        
        # Convert gext (S/cm^2) to an absolute conductance weight in microsiemens (uS)
        # The final scaling factor may need to be adjusted to match the old model's firing rate.
        # The term `rate * tau * 1e-3` approximates the average conductance.
        # Let's start with a value derived from your gext parameter.
        g_Abs_uS = config.current_state_params['gext_TC'] * area_cm2 * 1e6

        # The weight needs to be tuned. A good heuristic is to scale the absolute conductance
        # by a factor related to the average number of simultaneous events.
        # Let's start simple: the kick value.
        kick = 1.0 # From your getPoissonGating call
        
        # A reasonable weight is a fraction of the total absolute conductance. Let's start with a small value.
        # The exact value needs to be tuned to match the behavior of your old model.
        # A good starting point for the NetCon weight (in uS) is often a small number like 0.005 to 0.05.
        weight_uS = g_Abs_uS # (uS) <--- THIS IS THE KEY PARAMETER TO TUNE
        
        self.nc.weight[0] = weight_uS

        # print(f"Attached NetStim-based Poisson noise to TC cell {self.gid}")

    def createSynapses(self):
        pass

    def createInDegVars(self):
        pass

class RECell(Cell):
    def defineGeometry(self):
        self.soma.L = self.soma.diam = np.sqrt(1.43e4/np.pi) # one-compartment of 1.43e4 um2
        self.soma.nseg = 1

    def defineBiophysics(self):
        self.soma.insert('RE')

    def applyPoisson(self):
        """Applies CVode-safe Poisson noise using NetStim and ExpSyn."""
        
        # 1. Set parameters for the noise from your config file
        rate = config.baseline_RE  # e.g., 40 Hz
        
        # 2. Create the NetStim spike generator
        # We attach it to `self` to prevent it from being garbage collected
        self.ns = h.NetStim()
        self.ns.interval = 1000.0 / rate  # Convert rate (Hz) to mean interval (ms)
        self.ns.number = 1e9                     # Large number for continuous firing
        self.ns.noise = 1                 # A value of 1 makes the intervals exponential (a Poisson process)
        self.ns.start = 0                 # Start generating spikes at t=0

        # 3. Create a simple exponential synapse to receive the noise
        # ExpSyn is a built-in NEURON mechanism that does exactly what your Python code did.
        self.syn = h.ExpSyn(self.soma(0.5))
        self.syn.e = config.Eext_RE    # Reversal potential for the synapse
        
        # This is the time constant of the conductance decay.
        # It should be the SAME as the `tau` you used in your python script.
        self.syn.tau = config.tau_RE   # e.g., 2 ms

        # 4. Connect the spike generator to the synapse with a NetCon
        self.nc = h.NetCon(self.ns, self.syn)
        self.nc.delay = 0.1  # Use a minimal, non-zero delay for numerical stability

        # The 'weight' of the NetCon corresponds to the peak conductance increase (g_max) 
        # for a single noisy event. This is the equivalent of the 'kick' in your python script.
        # We will calculate a good starting value from your old `gext_PYdr` parameter.

        area_um2 = self.soma.L * np.pi * self.soma.diam
        area_cm2 = area_um2 * 1e-8
        
        # Convert gext (S/cm^2) to an absolute conductance weight in microsiemens (uS)
        # The final scaling factor may need to be adjusted to match the old model's firing rate.
        # The term `rate * tau * 1e-3` approximates the average conductance.
        # Let's start with a value derived from your gext parameter.
        g_Abs_uS = config.current_state_params['gext_RE'] * area_cm2 * 1e6

        # The weight needs to be tuned. A good heuristic is to scale the absolute conductance
        # by a factor related to the average number of simultaneous events.
        # Let's start simple: the kick value.
        kick = 1.0 # From your getPoissonGating call
        
        # A reasonable weight is a fraction of the total absolute conductance. Let's start with a small value.
        # The exact value needs to be tuned to match the behavior of your old model.
        # A good starting point for the NetCon weight (in uS) is often a small number like 0.005 to 0.05.
        weight_uS = g_Abs_uS # (uS) <--- THIS IS THE KEY PARAMETER TO TUNE
        
        self.nc.weight[0] = weight_uS

        # print(f"Attached NetStim-based Poisson noise to RE cell {self.gid}")

    def createSynapses(self):
        pass

    def createInDegVars(self):
        pass