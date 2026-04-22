'''
Defines network connectivity
'''
# --- import system libraries ---
import os
import sys
from datetime import datetime

# --- import the required libraries ---
import config
import pickle
from neuron import h
import numpy as np
from netcon import netcon_nearest_neighbors
# import matplotlib.pyplot as pyplot
h('load_file("stdgui.hoc")') #need this instead of import gui to get the simulation to be reproducible and not give an LFP flatline
from cells import Cell, PYsoCell, PYdrCell, INCell, TCCell, RECell   

class Net:
    """Creates network with prescribed number of each species of neurons (using parallelContext).
    Also ncludes methods to gather and plot spikes
    """
    def __init__(self, Npyso, Npydr, Ninh, Ntc, Nre):

        self.Npyso=Npyso                 #number of pyramidal cells
        self.Npydr=Npydr                 #number of pyramidal cellss
        self.Ninh=Ninh                   #number of inhibitory cortical cells
        self.Ntc=Ntc                     #number of thalamic cells
        self.Nre=Nre                     #number of reticular cells

        self.N = Npyso+Npydr+Ninh+Ntc+Nre     # total number of cells in network
        self.cells = []                       # List of Cell objects in the net
        self.nclist = []                      # List of NetCon in the net
        self.tVec = h.Vector()                # spike time of all cells on this processor
        self.idVec = h.Vector()               # cell ids of spike times

        # # New: Add EEG recording vectors
        # self.eeg_time_vec = None
        # self.eeg_signal_vec = None
        # self.pydr_cells = []  # Store PYdr cell objects directly

        self.createNet()                      # Actually build the net
        
    def __del__(self):
        config.pc.gid_clear()
        
    def createNet(self):
        """Create, layout, and connect N cells."""
        self.setGids() #### set global ids (gids), used to connect cells
        # self.createCells()

        try:
            self.createCells()
            print(f"Process {config.idhost}: Finished creating cells")
        except Exception as e:
            print(f"Process {config.idhost}: Error in createCells: {e}")
            raise

        # ---> Call the new sanity check function here <---
        self.check_network_consistency()

        config.pc.barrier()

        self.connectCells()
        self.createPoissonInput() 
        self.createStims() 
        self.createIClamps()

        # setup parallel transfers for continuous pointers
        config.pc.setup_transfer()

    def check_network_consistency(self):
        """
        The ultimate diagnostic. Each rank reports the GID, type, and exact
        list of mechanisms for each of its local cells. Rank 0 verifies consistency.
        """
        local_cells_info = []
        for cell in self.cells:
            cell_type = type(cell).__name__
            
            # --- Inspect the soma to get the TRUE list of mechanisms ---
            mechanisms = set()
            if hasattr(cell, 'soma'):
                for seg in cell.soma:
                    for mech in seg:
                        mechanisms.add(mech.name())
            
            local_cells_info.append({'gid': cell.gid, 'type': cell_type, 'mechs': sorted(list(mechanisms))})

        # Gather all lists onto rank 0
        all_ranks_info = config.pc.py_gather(local_cells_info, 0)
        
        if config.idhost == 0:
            print("\n--- RUNNING ULTIMATE NETWORK CONSISTENCY CHECK ---", flush=True)
            
            global_cell_map = {}
            found_error = False

            for rank_info in all_ranks_info:
                for cell_info in rank_info:
                    gid = cell_info['gid']
                    if gid in global_cell_map:
                        print(f"!!! FATAL ERROR: Duplicate GID {gid} detected!", flush=True)
                        found_error = True
                    global_cell_map[gid] = cell_info

            # --- Perform the detailed check ---
            base_mechs = {} # Store the "correct" set of mechanisms for each cell type
            for gid in sorted(global_cell_map.keys()):
                info = global_cell_map[gid]
                cell_type = info['type']
                mechs = str(info['mechs']) # Convert list to string for easy comparison

                if cell_type not in base_mechs:
                    # This is the first time we've seen this cell type. Assume it's the correct template.
                    base_mechs[cell_type] = mechs
                    print(f"Template for {cell_type}: GID {gid}, Mechs: {mechs}", flush=True)
                else:
                    # Compare this cell's mechanisms to the template for its type
                    if mechs != base_mechs[cell_type]:
                        print(f"!!! INCONSISTENCY FOUND !!!", flush=True)
                        print(f"  > GID {gid} ({cell_type}) has mechanisms: {mechs}", flush=True)
                        print(f"  > Expected mechanisms based on first instance: {base_mechs[cell_type]}", flush=True)
                        found_error = True

            if not found_error:
                print("--- ULTIMATE NETWORK CONSISTENCY CHECK PASSED! ---\n", flush=True)
            else:
                print("--- ULTIMATE NETWORK CONSISTENCY CHECK FAILED! ---\n", flush=True)
                config.pc.abort(-1)
                
        sys.stdout.flush()
        config.pc.barrier()
        
    def setGids(self):
        self.gidList = []
        self.pyso_gidList = []
        self.pydr_gidList = []
        self.inh_gidList = []
        self.tc_gidList = []
        self.re_gidList = []

        #### Round-robin counting. Each host as an id from 0 to nhost - 1.
        for i in range(config.idhost, self.N, config.nhost): #start with idhost, count by nhost until you get to total number of neurons (so if idhost=2 and nhost=4, the list contains [2,6,10...])
            self.gidList.append(i)
            if i<self.Npyso:
                self.pyso_gidList.append(i)
            elif(self.Npyso <= i < self.Npyso+self.Npydr):
                self.pydr_gidList.append(i)
            elif(self.Npyso+self.Npydr <= i < self.Npyso+self.Npydr+self.Ninh):
                self.inh_gidList.append(i)
            elif(self.Npyso+self.Npydr+self.Ninh <= i < self.Npyso+self.Npydr+self.Ninh+self.Ntc):
                self.tc_gidList.append(i)
            else:
                self.re_gidList.append(i)
 
    def createCells(self):
        """Create and layout cells (in this host) in the network."""
        self.cells = []
        
        for gid in self.gidList:  # Loop over cells in this node/host
            cell = None  # Initialize to None
            
            if gid < self.Npyso:
                cell = PYsoCell(gid)  # create pyramidal cell soma
                
            elif self.Npyso <= gid < self.Npyso + self.Npydr:
                cell = PYdrCell(gid)  # create pyramidal cell dendrite
                
            elif self.Npyso + self.Npydr <= gid < self.Npyso + self.Npydr + self.Ninh:
                cell = INCell(gid)  # create inhibitory cell
                
            elif self.Npyso + self.Npydr + self.Ninh <= gid < self.Npyso + self.Npydr + self.Ninh + self.Ntc:
                # Create TC cells only if thalamus is included in this condition
                if config.current_state_params['include_thalamus']:
                    cell = TCCell(gid)
                else:
                    print(f'Skipping TC cell {gid} (thalamus not included in {config.CONDITION} condition)')
                    continue  # Skip to next gid
                    
            else:  # RE cells
                # Create RE cells only if thalamus is included in this condition
                if config.current_state_params['include_thalamus']:
                    cell = RECell(gid)
                else:
                    print(f'Skipping RE cell {gid} (thalamus not included in {config.CONDITION} condition)')
                    continue  # Skip to next gid
            
            # Only proceed if a cell was created
            if cell is not None:
                self.cells.append(cell)  # add cell object to net cell list
                cell.associateGid()  # associate gid to each cell
                config.pc.spike_record(cell.gid, self.tVec, self.idVec)  # Record spikes

                # Setup voltage recording ONLY if the flag is passed
                if config.record_voltages:
                    cell.setRecording()
                
                # print(f'Created cell {gid} ({type(cell).__name__}) on host {config.idhost} out of {config.nhost}')
                # print(f'config.pc.gid2cell({gid}): {config.pc.gid2cell(gid)}')

    def connectCells(self):
        """Connect cells based on current condition. Note that this method assumes that there Npydr=100, Npyso=100, Ninh=20,
        Ntc=20, and Nre=20. It may not work correctly for other values."""

        if config.idhost == 0:
            print(f"Connecting cells for condition: {config.CONDITION}")
            print(f"Include thalamus: {config.current_state_params['include_thalamus']}")
            print(f"Include thalamocortical connections: {config.current_state_params['include_thalamocortical_connections']}")

        ### Setup the parallel transfer of presynaptic voltage on different hosts to synapses present on other hosts ###
        ### For more details refer to: https://nrn.readthedocs.io/en/latest/python/modelspec/programmatic/network/parcon.html#parallel-transfer ###

        # # Set up the transfer of presynaptic voltages for each cell type using gid_offset to distinguish it from the regular cell gids.    
        gid_offset = self.N

        # Transfer the PYso presynaptic voltages
        for pyso_gid in self.pyso_gidList:
            cell = config.pc.gid2cell(pyso_gid)
            config.pc.source_var(cell.soma(0.5)._ref_v, pyso_gid + gid_offset, sec=cell.soma)

        # Transfer the PYdr presynaptic voltages
        for pydr_gid in self.pydr_gidList:
            cell = config.pc.gid2cell(pydr_gid)
            config.pc.source_var(cell.soma(0.5)._ref_v, pydr_gid + gid_offset, sec=cell.soma)

        # Transfer the IN presynaptic voltages
        for inh_gid in self.inh_gidList:
            cell = config.pc.gid2cell(inh_gid)
            config.pc.source_var(cell.soma(0.5)._ref_v, inh_gid + gid_offset, sec=cell.soma)

        if config.current_state_params['include_thalamus']:
            # Transfer the TC presynaptic voltages
            for tc_gid in self.tc_gidList:
                cell = config.pc.gid2cell(tc_gid)
                config.pc.source_var(cell.soma(0.5)._ref_v, tc_gid + gid_offset, sec=cell.soma)

            # Transfer the RE presynaptic voltages
            for re_gid in self.re_gidList:
                cell = config.pc.gid2cell(re_gid)
                config.pc.source_var(cell.soma(0.5)._ref_v, re_gid + gid_offset, sec=cell.soma)

        # ================================================================
        # CORTICAL CONNECTIONS (Always included)
        # ================================================================

        ### More robust implementation of compartmental connections between PYso and PYdr ###
        if config.idhost == 0:
            print("Establishing PYso -> PYdr and PYdr -> PYso direct compartmental connections ...")

        # --- PYso -> PYdr direct connections via ICOM ---

        PYso_PYdr_conn = np.eye(self.Npyso)
        # Extract sources and targets from the connectivity matrix
        PYso_sources, PYdr_targets= np.where(PYso_PYdr_conn == 1)  # Get indices where conn[i, j] = 1
        self.PYso_PYdr_ICOM_synapses = []
        for src, tgt in zip(PYso_sources, PYdr_targets):
            pyso_src_gid = src
            pydr_tgt_gid = tgt + self.Npyso
            if config.pc.gid_exists(pydr_tgt_gid):
                pydr_cell = config.pc.gid2cell(pydr_tgt_gid)
                icom_syn = h.ICOM(pydr_cell.soma(0.5))

                # Set synapse parameters
                icom_syn.area_post = pydr_cell.soma(0.5).area() * 1e-8 # um2 to cm2

                # Store the synapses in the list
                self.PYso_PYdr_ICOM_synapses.append(icom_syn)

                # Presynaptic cell is local or remote
                config.pc.target_var(icom_syn, icom_syn._ref_v_pre, pyso_src_gid + gid_offset)
                # print(f"Connected PYso {pyso_src_gid} to PYdr {pydr_tgt_gid} via ICOM")

        # --- PYdr -> PYso direct connections via ICOM and IKNa ---

        PYdr_PYso_conn = np.eye(self.Npydr)
        # Extract sources and targets from the connectivity matrix
        PYdr_sources, PYso_targets= np.where(PYdr_PYso_conn == 1)  # Get indices where conn[i, j] = 1
        self.PYdr_PYso_ICOM_synapses = []
        self.PYdr_PYso_IKNa_synapses = []
        for src, tgt in zip(PYdr_sources, PYso_targets):
            pydr_src_gid = src + self.Npyso
            pyso_tgt_gid = tgt
            if config.pc.gid_exists(pyso_tgt_gid):
                pyso_cell = config.pc.gid2cell(pyso_tgt_gid)
                icom_syn = h.ICOM(pyso_cell.soma(0.5))
                ikna_syn = h.IKNa(pyso_cell.soma(0.5))

                # Set synapse parameters
                icom_syn.area_post = pyso_cell.soma(0.5).area() * 1e-8 # um2 to cm2
                ikna_syn.area_post = pyso_cell.soma(0.5).area() * 1e-8 # um2 to cm2
                ikna_syn.gKNa = config.current_state_params['gKNa'] # S/cm^2

                # Store the synapses in the list
                self.PYdr_PYso_ICOM_synapses.append(icom_syn)
                self.PYdr_PYso_IKNa_synapses.append(ikna_syn)

                # Presynaptic cell is local or remote
                config.pc.target_var(icom_syn, icom_syn._ref_v_pre, pydr_src_gid + gid_offset)
                config.pc.target_var(ikna_syn, ikna_syn._ref_v_pre, pydr_src_gid + gid_offset)
                # print(f"Connected PYdr {pydr_src_gid} to PYso {pyso_tgt_gid} via ICOM and IKNa")

        ### AMPAergic, NMDAergic, and GABAAergic intracortical connections ###

        Radius = 10
        # We want to remove recurrent connections for PYso -> PYdr and PYdr -> PYso to avoid self-excitation and self-inhibition.
        Remove_Recurrent_Bool = True 
        if config.idhost == 0:
            print("Establishing intracortical connections ...")

        # 1) Connect PYso -> PYdr with AMPAergic synapses
        conn1 = netcon_nearest_neighbors(2 * Radius, self.Npyso, self.Npydr, Remove_Recurrent_Bool)
        sources1, targets1 = np.where(conn1 == 1)  # Get indices where conn[i, j] = 1
        self.AMPA_PYso_PYdr_synapses = []
        self.AMPA_PYso_PYdr_netcons = []
        self.AMPA_PYso_PYdr_per_PYdr = [[] for _ in range(self.Npydr)]
        for src, tgt in zip(sources1, targets1):
            pyso_src_gid = src
            pydr_tgt_gid = tgt + self.Npyso
            if config.pc.gid_exists(pydr_tgt_gid):
                pydr_cell = config.pc.gid2cell(pydr_tgt_gid)
                syn = h.AMPAD(pydr_cell.soma(0.5))
                syn.Npre = self.Npyso
                syn.Npost = self.Npydr
                syn.gAMPA = config.current_state_params['gAMPA_PYso_PYdr']
                syn.EAMPA = 0
                syn.alpha = 3.48
                syn.tau = 2
                syn.deprFactor = config.deprFactor
                syn.tauRes = config.tauRes
                syn.release_time = config.release_time
                syn.radius = 10
                syn.remove_recurrent_bool = 1 if Remove_Recurrent_Bool else 0
                syn.area_post = pydr_cell.soma(0.5).area() * 1e-8
                self.AMPA_PYso_PYdr_synapses.append(syn)
                self.AMPA_PYso_PYdr_per_PYdr[tgt].append(syn)

                # ==========================================================
                # CRITICAL: Use OFFSET GID for voltage transfer (continuous)
                # ==========================================================
                
                # 1. The Analog Link (For Graded Activation - sAMPA/xNMDA)
                # This drives the continuous Equation 1 in DERIVATIVE block
                config.pc.target_var(syn, syn._ref_v_pre, pyso_src_gid + gid_offset)

                # 2. The Digital Link (For Depression - res_AMPA)
                # This drives the NET_RECEIVE block to multiply by deprFactor
                # Note: We MUST create a NetCon.
                nc = config.pc.gid_connect(pyso_src_gid, syn)
                nc.threshold = config.threshold
                nc.delay = 0.1      # Set a minimal delay to avoid zero-delay issues for parallel simulation
                nc.weight[0] = 0     # Weight is unused in our MOD, so 0 is fine

                self.AMPA_PYso_PYdr_netcons.append(nc)
                # print(f"Process {config.idhost}: Connected PYso {pyso_src_gid} to PYdr {pydr_tgt_gid} with AMPA")

        # 2) Connect PYso -> PYdr with NMDAergic synapses
        conn2 = netcon_nearest_neighbors(2 * Radius, self.Npyso, self.Npydr, Remove_Recurrent_Bool)
        sources2, targets2 = np.where(conn2 == 1)  # Get indices where conn[i, j] = 1
        self.NMDA_PYso_PYdr_synapses = []
        self.NMDA_PYso_PYdr_netcons = []
        self.NMDA_PYso_PYdr_per_PYdr = [[] for _ in range(self.Npydr)]
        for src, tgt in zip(sources2, targets2):
            pyso_src_gid = src
            pydr_tgt_gid = tgt + self.Npyso
            if config.pc.gid_exists(pydr_tgt_gid):
                pydr_cell = config.pc.gid2cell(pydr_tgt_gid)
                syn = h.NMDAD(pydr_cell.soma(0.5))
                syn.Npre = self.Npyso
                syn.Npost = self.Npydr
                syn.gNMDA = config.current_state_params['gNMDA_PYso_PYdr']
                syn.ENMDA = 0
                syn.alphaS = 0.5
                syn.tauS = 100
                syn.alphaX = 3.48
                syn.tauX = 2
                syn.deprFactor = config.deprFactor
                syn.tauRes = config.tauRes
                syn.release_time = config.release_time
                syn.radius = 10
                syn.remove_recurrent_bool = 1 if Remove_Recurrent_Bool else 0
                syn.area_post = pydr_cell.soma(0.5).area() * 1e-8
                self.NMDA_PYso_PYdr_synapses.append(syn)
                self.NMDA_PYso_PYdr_per_PYdr[tgt].append(syn)

                # ==========================================================
                # CRITICAL: Use OFFSET GID for voltage transfer (continuous)
                # ==========================================================
                
                # 1. The Analog Link (For Graded Activation - sAMPA/xNMDA)
                # This drives the continuous Equation 1 in DERIVATIVE block
                config.pc.target_var(syn, syn._ref_v_pre, pyso_src_gid + gid_offset)

                # 2. The Digital Link (For Depression - res_AMPA)
                # This drives the NET_RECEIVE block to multiply by deprFactor
                # Note: We MUST create a NetCon.
                nc = config.pc.gid_connect(pyso_src_gid, syn)
                nc.threshold = config.threshold
                nc.delay = 0.1      # Set a minimal delay to avoid zero-delay issues for parallel simulation
                nc.weight[0] = 0     # Weight is unused in our MOD, so 0 is fine

                self.NMDA_PYso_PYdr_netcons.append(nc)
                # print(f"Process {config.idhost}: Connected PYso {pyso_src_gid} to PYdr {pydr_tgt_gid} with NMDA")

        # --- no need to remove recurrent connections for PYso -> IN and IN -> PYso since they are different cell populations ---
        Remove_Recurrent_Bool = False

        # 3) Connect PYso -> IN with AMPAergic synapses
        conn3 = netcon_nearest_neighbors(2 * Radius, self.Npyso, self.Ninh, Remove_Recurrent_Bool)
        sources3, targets3 = np.where(conn3 == 1)  # Get indices where conn[i, j] = 1
        self.AMPA_PYso_IN_synapses = []
        self.AMPA_PYso_IN_netcons = []
        self.AMPA_PYso_IN_per_IN = [[] for _ in range(self.Ninh)]
        for src, tgt in zip(sources3, targets3):
            pyso_src_gid = src
            inh_tgt_gid = tgt + self.Npyso + self.Npydr
            if config.pc.gid_exists(inh_tgt_gid):
                inh_cell = config.pc.gid2cell(inh_tgt_gid)
                syn = h.AMPAD(inh_cell.soma(0.5))
                syn.Npre = self.Npyso
                syn.Npost = self.Ninh
                syn.gAMPA = config.current_state_params['gAMPA_PYso_IN']
                syn.EAMPA = 0
                syn.alpha = 3.48
                syn.tau = 2
                syn.deprFactor = config.deprFactor
                syn.tauRes = config.tauRes
                syn.release_time = config.release_time
                syn.radius = 10
                syn.remove_recurrent_bool = 1 if Remove_Recurrent_Bool else 0
                syn.area_post = inh_cell.soma(0.5).area() * 1e-8
                self.AMPA_PYso_IN_synapses.append(syn)
                self.AMPA_PYso_IN_per_IN[tgt].append(syn)

                # ==========================================================
                # CRITICAL: Use OFFSET GID for voltage transfer (continuous)
                # ==========================================================
                
                # 1. The Analog Link (For Graded Activation - sAMPA/xNMDA)
                # This drives the continuous Equation 1 in DERIVATIVE block
                config.pc.target_var(syn, syn._ref_v_pre, pyso_src_gid + gid_offset)

                # 2. The Digital Link (For Depression - res_AMPA)
                # This drives the NET_RECEIVE block to multiply by deprFactor
                # Note: We MUST create a NetCon.
                nc = config.pc.gid_connect(pyso_src_gid, syn)
                nc.threshold = config.threshold
                nc.delay = 0.1      # Set a minimal delay to avoid zero-delay issues for parallel simulation
                nc.weight[0] = 0     # Weight is unused in our MOD, so 0 is fine

                self.AMPA_PYso_IN_netcons.append(nc)
                # print(f"Process {config.idhost}: Connected PYso {pyso_src_gid} to IN {inh_tgt_gid} with AMPA")

        # 4) Connect PYso -> IN with NMDAergic synapses
        conn4 = netcon_nearest_neighbors(2 * Radius, self.Npyso, self.Ninh, Remove_Recurrent_Bool)
        sources4, targets4 = np.where(conn4 == 1)  # Get indices where conn[i, j] = 1
        self.NMDA_PYso_IN_synapses = []
        self.NMDA_PYso_IN_netcons = []
        self.NMDA_PYso_IN_per_IN = [[] for _ in range(self.Ninh)]
        for src, tgt in zip(sources4, targets4):
            pyso_src_gid = src
            inh_tgt_gid = tgt + self.Npyso + self.Npydr
            if config.pc.gid_exists(inh_tgt_gid):
                inh_cell = config.pc.gid2cell(inh_tgt_gid)
                syn = h.NMDAD(inh_cell.soma(0.5))
                syn.Npre = self.Npyso
                syn.Npost = self.Ninh
                syn.gNMDA = config.current_state_params['gNMDA_PYso_IN']
                syn.ENMDA = 0
                syn.alphaS = 0.5
                syn.tauS = 100
                syn.alphaX = 3.48
                syn.tauX = 2
                syn.deprFactor = config.deprFactor
                syn.tauRes = config.tauRes
                syn.release_time = config.release_time
                syn.radius = 10
                syn.remove_recurrent_bool = 1 if Remove_Recurrent_Bool else 0
                syn.area_post = inh_cell.soma(0.5).area() * 1e-8
                self.NMDA_PYso_IN_synapses.append(syn)
                self.NMDA_PYso_IN_per_IN[tgt].append(syn)

                # ==========================================================
                # CRITICAL: Use OFFSET GID for voltage transfer (continuous)
                # ==========================================================
                
                # 1. The Analog Link (For Graded Activation - sAMPA/xNMDA)
                # This drives the continuous Equation 1 in DERIVATIVE block
                config.pc.target_var(syn, syn._ref_v_pre, pyso_src_gid + gid_offset)

                # 2. The Digital Link (For Depression - res_AMPA)
                # This drives the NET_RECEIVE block to multiply by deprFactor
                # Note: We MUST create a NetCon.
                nc = config.pc.gid_connect(pyso_src_gid, syn)
                nc.threshold = config.threshold
                nc.delay = 0.1      # Set a minimal delay to avoid zero-delay issues for parallel simulation
                nc.weight[0] = 0     # Weight is unused in our MOD, so 0 is fine

                self.NMDA_PYso_IN_netcons.append(nc)
                # print(f"Process {config.idhost}: Connected PYso {pyso_src_gid} to IN {inh_tgt_gid} with NMDA")

        # 5) Connect IN -> PYso with GABAAergic synapses
        conn5 = netcon_nearest_neighbors(2 * Radius, self.Ninh, self.Npyso, Remove_Recurrent_Bool)
        sources5, targets5 = np.where(conn5 == 1)  # Get indices where conn[i, j] = 1
        self.GABAA_IN_PYso_synapses = []
        self.GABAA_IN_PYso_netcons = []
        self.GABAA_IN_PYso_per_PYso = [[] for _ in range(self.Npyso)]
        for src, tgt in zip(sources5, targets5):
            inh_src_gid = src + self.Npyso + self.Npydr
            pyso_tgt_gid = tgt
            if config.pc.gid_exists(pyso_tgt_gid):
                pyso_cell = config.pc.gid2cell(pyso_tgt_gid)
                syn = h.GABAAD(pyso_cell.soma(0.5))
                syn.Npre = self.Ninh
                syn.Npost = self.Npyso
                syn.gGABAA = config.current_state_params['gGABAA_IN_PYso']
                syn.EGABAA = -70
                syn.alpha = 1
                syn.tauGABAA = 5
                syn.deprFactor = config.deprFactor
                syn.tauRes = config.tauRes
                syn.release_time = config.release_time
                syn.radius = 10
                syn.remove_recurrent_bool = 1 if Remove_Recurrent_Bool else 0
                syn.propoCondMult = config.current_state_params['propoCondMult']
                syn.propoTauMult = config.current_state_params['propoTauMult']
                syn.area_post = pyso_cell.soma(0.5).area() * 1e-8
                self.GABAA_IN_PYso_synapses.append(syn)
                self.GABAA_IN_PYso_per_PYso[tgt].append(syn)

                # ==========================================================
                # CRITICAL: Use OFFSET GID for voltage transfer (continuous)
                # ==========================================================
                
                # 1. The Analog Link (For Graded Activation - sAMPA/xNMDA)
                # This drives the continuous Equation 1 in DERIVATIVE block
                config.pc.target_var(syn, syn._ref_v_pre, inh_src_gid + gid_offset)

                # 2. The Digital Link (For Depression - res_AMPA)
                # This drives the NET_RECEIVE block to multiply by deprFactor
                # Note: We MUST create a NetCon.
                nc = config.pc.gid_connect(inh_src_gid, syn)
                nc.threshold = config.threshold
                nc.delay = 0.1      # Set a minimal delay to avoid zero-delay issues for parallel simulation
                nc.weight[0] = 0     # Weight is unused in our MOD, so 0 is fine

                self.GABAA_IN_PYso_netcons.append(nc)
                # print(f"Process {config.idhost}: Connected IN {inh_src_gid} to PYso {pyso_tgt_gid} with GABAA")

        # --- We need to remove recurrent connections for IN -> IN to avoid self-inhibition ---
        Remove_Recurrent_Bool = True

        # 6) Connect IN -> IN with GABAAergic synapses
        conn6 = netcon_nearest_neighbors(2 * Radius, self.Ninh, self.Ninh, Remove_Recurrent_Bool)
        sources6, targets6 = np.where(conn6 == 1)  # Get indices where conn[i, j] = 1
        self.GABAA_IN_IN_synapses = []
        self.GABAA_IN_IN_netcons = []
        self.GABAA_IN_IN_per_IN = [[] for _ in range(self.Ninh)]
        for src, tgt in zip(sources6, targets6):
            inh_src_gid = src + self.Npyso + self.Npydr
            inh_tgt_gid = tgt + self.Npyso + self.Npydr
            if config.pc.gid_exists(inh_tgt_gid):
                inh_cell = config.pc.gid2cell(inh_tgt_gid)
                syn = h.GABAAD(inh_cell.soma(0.5))
                syn.Npre = self.Ninh
                syn.Npost = self.Ninh
                syn.gGABAA = 0.000825e-3
                syn.EGABAA = -70
                syn.alpha = 1
                syn.tauGABAA = 5
                syn.deprFactor = config.deprFactor
                syn.tauRes = config.tauRes
                syn.release_time = config.release_time
                syn.radius = 10
                syn.remove_recurrent_bool = 1 if Remove_Recurrent_Bool else 0
                syn.propoCondMult = config.current_state_params['propoCondMult']
                syn.propoTauMult = config.current_state_params['propoTauMult']
                syn.area_post = inh_cell.soma(0.5).area() * 1e-8
                self.GABAA_IN_IN_synapses.append(syn)
                self.GABAA_IN_IN_per_IN[tgt].append(syn)

                # ==========================================================
                # CRITICAL: Use OFFSET GID for voltage transfer (continuous)
                # ==========================================================
                
                # 1. The Analog Link (For Graded Activation - sAMPA/xNMDA)
                # This drives the continuous Equation 1 in DERIVATIVE block
                config.pc.target_var(syn, syn._ref_v_pre, inh_src_gid + gid_offset)

                # 2. The Digital Link (For Depression - res_AMPA)
                # This drives the NET_RECEIVE block to multiply by deprFactor
                # Note: We MUST create a NetCon.
                nc = config.pc.gid_connect(inh_src_gid, syn)
                nc.threshold = config.threshold
                nc.delay = 0.1      # Set a minimal delay to avoid zero-delay issues for parallel simulation
                nc.weight[0] = 0     # Weight is unused in our MOD, so 0 is fine

                self.GABAA_IN_IN_netcons.append(nc)
                # print(f"Process {config.idhost}: Connected IN {inh_src_gid} to IN {inh_tgt_gid} with GABAA")

        if config.current_state_params['include_thalamus']:

            ### Set up the Thalamic connections ###
            if config.idhost == 0:
                print("Establishing thalamic connections ...")

            # 7) Connect TC -> RE with AMPAergic synapses
            conn7 = np.ones((self.Ntc, self.Nre)) # all to all connections
            sources7, targets7 = np.where(conn7 == 1)  # Get indices where conn[i, j] = 1
            self.AMPA_TC_RE_synapses = []
            self.AMPA_TC_RE_per_RE = [[] for _ in range(self.Nre)]
            for src, tgt in zip(sources7, targets7):
                tc_src_gid = src + self.Npyso + self.Npydr + self.Ninh
                re_tgt_gid = tgt + self.Npyso + self.Npydr + self.Ninh + self.Ntc
                if config.pc.gid_exists(re_tgt_gid):
                    re_cell = config.pc.gid2cell(re_tgt_gid)
                    syn = h.AMPA(re_cell.soma(0.5))
                    syn.normalizing_factor = self.Ntc
                    syn.gAMPA = 0.4e-3  # S/cm^2
                    syn.EAMPA = 1 # mV
                    syn.alpha = 5 # ms^-1
                    syn.tauAMPA = 2 # ms
                    syn.area_post = re_cell.soma(0.5).area() * 1e-8 # um2 to cm2
                    self.AMPA_TC_RE_synapses.append(syn)
                    self.AMPA_TC_RE_per_RE[tgt].append(syn)
                    config.pc.target_var(syn, syn._ref_v_pre, tc_src_gid + gid_offset)
                    # print(f"Process {config.idhost}: Connected TC {tc_src_gid} to RE {re_tgt_gid} with AMPA")

            # 8) Connect RE -> RE with GABAAergic synapses
            conn8 = np.ones((self.Nre, self.Nre)) # all to all connections
            sources8, targets8 = np.where(conn8 == 1)  # Get indices where conn[i, j] = 1
            self.GABAA_RE_RE_synapses = []
            self.GABAA_RE_RE_per_RE = [[] for _ in range(self.Nre)]
            for src, tgt in zip(sources8, targets8):
                re_src_gid = src + self.Npyso + self.Npydr + self.Ninh + self.Ntc
                re_tgt_gid = tgt + self.Npyso + self.Npydr + self.Ninh + self.Ntc
                if config.pc.gid_exists(re_tgt_gid):
                    re_cell = config.pc.gid2cell(re_tgt_gid)
                    syn = h.GABAA(re_cell.soma(0.5))
                    syn.normalizing_factor = self.Nre
                    syn.gGABAA = 0.1e-3  # S/cm^2
                    syn.EGABAA = -80 # mV
                    syn.alpha = 2 # ms^-1
                    syn.tauGABAA = 5 # ms
                    syn.area_post = re_cell.soma(0.5).area() * 1e-8 # um2 to cm2
                    syn.propoCondMult = config.current_state_params['propoCondMult']
                    syn.propoTauMult = config.current_state_params['propoTauMult']
                    self.GABAA_RE_RE_synapses.append(syn)
                    self.GABAA_RE_RE_per_RE[tgt].append(syn)
                    config.pc.target_var(syn, syn._ref_v_pre, re_src_gid + gid_offset)
                    # print(f"Process {config.idhost}: Connected RE {re_src_gid} to RE {re_tgt_gid} with GABAA")

            # 9) Connect RE -> TC with GABAAergic synapses
            conn9 = np.ones((self.Nre, self.Ntc)) # all to all connections
            sources9, targets9 = np.where(conn9 == 1)  # Get indices where conn[i, j] = 1
            self.GABAA_RE_TC_synapses = []
            self.GABAA_RE_TC_per_TC = [[] for _ in range(self.Ntc)]
            for src, tgt in zip(sources9, targets9):
                re_src_gid = src + self.Npyso + self.Npydr + self.Ninh + self.Ntc
                tc_tgt_gid = tgt + self.Npyso + self.Npydr + self.Ninh
                if config.pc.gid_exists(tc_tgt_gid):
                    tc_cell = config.pc.gid2cell(tc_tgt_gid)
                    syn = h.GABAA(tc_cell.soma(0.5))
                    syn.normalizing_factor = self.Nre
                    syn.gGABAA = 0.1e-3  # S/cm^2
                    syn.EGABAA = -80 # mV
                    syn.alpha = 2 # ms^-1
                    syn.tauGABAA = 5 # ms
                    syn.area_post = tc_cell.soma(0.5).area() * 1e-8 # um2 to cm2
                    syn.propoCondMult = config.current_state_params['propoCondMult']
                    syn.propoTauMult = config.current_state_params['propoTauMult']
                    self.GABAA_RE_TC_synapses.append(syn)
                    self.GABAA_RE_TC_per_TC[tgt].append(syn)
                    config.pc.target_var(syn, syn._ref_v_pre, re_src_gid + gid_offset)
                    # print(f"Process {config.idhost}: Connected RE {re_src_gid} to TC {tc_tgt_gid} with GABAA")

            # 10) Connect RE -> TC with GABABergic synapses
            conn10 = np.ones((self.Nre, self.Ntc)) # all to all connections
            sources10, targets10 = np.where(conn10 == 1)  # Get indices where conn[i, j] = 1
            self.GABAB_RE_TC_synapses = []
            self.GABAB_RE_TC_per_TC = [[] for _ in range(self.Ntc)]
            for src, tgt in zip(sources10, targets10):
                re_src_gid = src + self.Npyso + self.Npydr + self.Ninh + self.Ntc
                tc_tgt_gid = tgt + self.Npyso + self.Npydr + self.Ninh
                if config.pc.gid_exists(tc_tgt_gid):
                    tc_cell = config.pc.gid2cell(tc_tgt_gid)
                    syn = h.GABAB(tc_cell.soma(0.5))
                    syn.normalizing_factor = self.Nre
                    syn.gGABAB = 0.001e-3  # S/cm^2
                    syn.EGABAB = -95 # mV
                    syn.K1 = 0.5 # ms-1
                    syn.K2 = 0.0012 # ms-1
                    syn.K3 = 0.18 # ms-1
                    syn.K4 = 0.034 # ms-1
                    syn.area_post = tc_cell.soma(0.5).area() * 1e-8 # um2 to cm2
                    self.GABAB_RE_TC_synapses.append(syn)
                    self.GABAB_RE_TC_per_TC[tgt].append(syn)
                    config.pc.target_var(syn, syn._ref_v_pre, re_src_gid + gid_offset)
                    # print(f"Process {config.idhost}: Connected RE {re_src_gid} to TC {tc_tgt_gid} with GABAB")

        if config.current_state_params['include_thalamocortical_connections']:

            ### Set up the thalamocortical connections ###
            Radius = 10
            Remove_Recurrent_Bool = False
            if config.idhost == 0:
                print("Establishing thalamocortical connections ...")

            # 11) Connect TC -> PYdr with AMPAergic synapses
            conn11 = netcon_nearest_neighbors(2 * Radius, self.Ntc, self.Npydr, Remove_Recurrent_Bool)
            sources11, targets11 = np.where(conn11 == 1)  # Get indices where conn[i, j] = 1
            self.AMPA_TC_PYdr_synapses = []
            self.AMPA_TC_PYdr_per_PYdr = [[] for _ in range(self.Npydr)]
            for src, tgt in zip(sources11, targets11):
                tc_src_gid = src + self.Npyso + self.Npydr + self.Ninh
                pydr_tgt_gid = tgt + self.Npyso
                if config.pc.gid_exists(pydr_tgt_gid):
                    pydr_cell = config.pc.gid2cell(pydr_tgt_gid)
                    syn = h.AMPA(pydr_cell.soma(0.5))
                    syn.normalizing_factor = np.clip((2 * Radius + (1 - int(Remove_Recurrent_Bool))) / (self.Npydr / self.Ntc), 0, self.Ntc)
                    syn.gAMPA = config.current_state_params['gAMPA_TC_PYdr']
                    syn.EAMPA = 1
                    syn.alpha = 5
                    syn.tauAMPA = 2
                    syn.area_post = pydr_cell.soma(0.5).area() * 1e-8
                    self.AMPA_TC_PYdr_synapses.append(syn)
                    self.AMPA_TC_PYdr_per_PYdr[tgt].append(syn)
                    config.pc.target_var(syn, syn._ref_v_pre, tc_src_gid + gid_offset)
                    # print(f"Process {config.idhost}: Connected TC {tc_src_gid} to PYdr {pydr_tgt_gid} with AMPA")

            # 12) Connect TC -> IN with AMPAergic synapses
            conn12 = netcon_nearest_neighbors(2 * Radius, self.Ntc, self.Ninh, Remove_Recurrent_Bool)
            sources12, targets12 = np.where(conn12 == 1)  # Get indices where conn[i, j] = 1
            self.AMPA_TC_IN_synapses = []
            self.AMPA_TC_IN_per_IN = [[] for _ in range(self.Ninh)]
            for src, tgt in zip(sources12, targets12):
                tc_src_gid = src + self.Npyso + self.Npydr + self.Ninh
                inh_tgt_gid = tgt + self.Npyso + self.Npydr
                if config.pc.gid_exists(inh_tgt_gid):
                    inh_cell = config.pc.gid2cell(inh_tgt_gid)
                    syn = h.AMPA(inh_cell.soma(0.5))
                    syn.normalizing_factor = np.clip((2 * Radius + (1 - int(Remove_Recurrent_Bool))) / (self.Ninh / self.Ntc), 0, self.Ntc)
                    syn.gAMPA = config.current_state_params['gAMPA_TC_IN']
                    syn.EAMPA = 1
                    syn.alpha = 5
                    syn.tauAMPA = 2
                    syn.area_post = inh_cell.soma(0.5).area() * 1e-8
                    self.AMPA_TC_IN_synapses.append(syn)
                    self.AMPA_TC_IN_per_IN[tgt].append(syn)
                    config.pc.target_var(syn, syn._ref_v_pre, tc_src_gid + gid_offset)
                    # print(f"Process {config.idhost}: Connected TC {tc_src_gid} to IN {inh_tgt_gid} with AMPA")

            # 13) Connect PYso -> TC with AMPAergic synapses
            conn13 = netcon_nearest_neighbors(2 * Radius, self.Npyso, self.Ntc, Remove_Recurrent_Bool)
            sources13, targets13 = np.where(conn13 == 1)  # Get indices where conn[i, j] = 1
            self.AMPA_PYso_TC_synapses = []
            self.AMPA_PYso_TC_per_TC = [[] for _ in range(self.Ntc)]
            for src, tgt in zip(sources13, targets13):
                pyso_src_gid = src
                tc_tgt_gid = tgt + self.Npyso + self.Npydr + self.Ninh
                if config.pc.gid_exists(tc_tgt_gid):
                    tc_cell = config.pc.gid2cell(tc_tgt_gid)
                    syn = h.AMPA(tc_cell.soma(0.5))
                    syn.normalizing_factor = np.clip((2 * Radius + (1 - int(Remove_Recurrent_Bool))) / (self.Ntc / self.Npyso), 0, self.Npyso)
                    syn.gAMPA = config.current_state_params['gAMPA_PYso_TC']
                    syn.EAMPA = 1
                    syn.alpha = 5
                    syn.tauAMPA = 2
                    syn.area_post = tc_cell.soma(0.5).area() * 1e-8
                    self.AMPA_PYso_TC_synapses.append(syn)
                    self.AMPA_PYso_TC_per_TC[tgt].append(syn)
                    config.pc.target_var(syn, syn._ref_v_pre, pyso_src_gid + gid_offset)
                    # print(f"Process {config.idhost}: Connected PYso {pyso_src_gid} to TC {tc_tgt_gid} with AMPA")

            # 14) Connect PYso -> RE with AMPAergic synapses
            conn14 = netcon_nearest_neighbors(2 * Radius, self.Npyso, self.Nre, Remove_Recurrent_Bool)
            sources14, targets14 = np.where(conn14 == 1)  # Get indices where conn[i, j] = 1
            self.AMPA_PYso_RE_synapses = []
            self.AMPA_PYso_RE_per_RE = [[] for _ in range(self.Nre)]
            for src, tgt in zip(sources14, targets14):
                pyso_src_gid = src
                re_tgt_gid = tgt + self.Npyso + self.Npydr + self.Ninh + self.Ntc
                if config.pc.gid_exists(re_tgt_gid):
                    re_cell = config.pc.gid2cell(re_tgt_gid)
                    syn = h.AMPA(re_cell.soma(0.5))
                    syn.normalizing_factor = np.clip((2 * Radius + (1 - int(Remove_Recurrent_Bool))) / (self.Nre / self.Npyso), 0, self.Npyso)
                    syn.gAMPA = config.current_state_params['gAMPA_PYso_RE']
                    syn.EAMPA = 1
                    syn.alpha = 5
                    syn.tauAMPA = 2
                    syn.area_post = re_cell.soma(0.5).area() * 1e-8
                    self.AMPA_PYso_RE_synapses.append(syn)
                    self.AMPA_PYso_RE_per_RE[tgt].append(syn)
                    config.pc.target_var(syn, syn._ref_v_pre, pyso_src_gid + gid_offset)
                    # print(f"Process {config.idhost}: Connected PYso {pyso_src_gid} to RE {re_tgt_gid} with AMPA")

    def createPoissonInput(self):
        """Set up background Poisson inputs/noise for PYdr, TC, and RE cells."""
        pass    
            
    def createStims(self):
        if config.stim_type == 'none':
            return

        if config.idhost == 0:    
            print(f"Applying {config.stim_type} stimulus starting at {config.stim_start} ms for {config.CONDITION} condition...")
        
        # Keep references to IClamps so they aren't garbage collected
        self.iclamp_list = []
        
        for pydr_gid in self.pydr_gidList:
            if config.pc.gid_exists(pydr_gid):
                cell = config.pc.gid2cell(pydr_gid)
                
                # 1. Create the current clamp
                stim = h.IClamp(cell.soma(0.5))
                stim.delay = config.stim_start
                stim.dur = config.stim_dur
                
                # 2. Calculate the surface area in cm^2
                area_cm2 = cell.soma(0.5).area() * 1e-8

                J_uA_cm2 = 0.0 # Default value, will be overwritten based on stim_type
                # 3. Determine the current density (J) in uA/cm2
                if config.stim_type == 'sync':
                    J_uA_cm2 = 1.0 
                elif config.stim_type == 'desync':
                    # Seed with gid to ensure parallel reproducibility
                    np.random.seed(int(pydr_gid + config.randSeed))
                    J_uA_cm2 = np.random.uniform(-1.0, 1.0)
                    
                # 4. Convert to nA: I(nA) = J(uA/cm2) * Area(cm2) * 1000
                amp_nA = J_uA_cm2 * area_cm2 * 1000.0
                stim.amp = amp_nA
                
                self.iclamp_list.append(stim)
    
    def createIClamps(self):
        pass
        
    def setCellLocations(self):
        '''set cortical cell locations for calculating distance from recording 
        electrode. Also set pointers necessary for xtra and extracellular 
        mechanisms to work. Note we are only interested in placing PYR and INH
        cells because RE and TC cells do not contribute to LFP.'''
        pass

    def place_cell(self):
        '''places zeroth 3D coordinate of soma of 'cell' at a specified x,y,z coordinate, maintaining default orientation of cell 
        (cells appear to be parallel to the x-axis) (This code is copied from HFO recode version 12)'''
        pass
    
    def gatherSpikes(self):
        """Gather spikes from all nodes/hosts"""
        if config.idhost==0: print('Gathering spikes ...')
        
        data = [None]*config.nhost
        data[0] = {'tVec': self.tVec, 'idVec': self.idVec}
        config.pc.barrier()
        gather=config.pc.py_alltoall(data)
        config.pc.barrier()
        self.tVecAll = [] 
        self.idVecAll = [] 
        if config.idhost==0:
            for d in gather:
                self.tVecAll.extend(list(d['tVec']))
                self.idVecAll.extend(list(d['idVec']))

    def gatherVoltages(self):
        """Gather voltage traces from all nodes to rank 0."""
        if config.idhost == 0: 
            print('Gathering voltages ...')
            
        # Extract local voltages into a dictionary mapping GID to numpy array
        local_v = {}
        for cell in self.cells:
            # Appending .copy() ensures safe serialization across MPI!
            local_v[cell.gid] = cell.soma_v_vec.as_numpy().copy() 
            
        # Get the global time vector, ensuring a safe copy
        local_t = self.cells[0].tVec.as_numpy().copy() if self.cells else np.array([])

        config.pc.barrier()
        # py_gather sends data from all ranks to rank 0
        gather_v = config.pc.py_gather(local_v, 0)
        gather_t = config.pc.py_gather(local_t, 0)
        config.pc.barrier()
        
        self.all_voltages = {}
        self.t_vec_global = np.array([])
        
        if config.idhost == 0:
            # Reconstruct the global dictionary of all cell voltages
            for d in gather_v:
                self.all_voltages.update(d)
                
            # Grab the time vector from the first rank that had cells
            for t in gather_t:
                if len(t) > 0:
                    self.t_vec_global = t
                    break
                
    def gatherLFP(self):
        '''Gather LFP waveforms from all nodes/hosts'''
        pass

    def gatherEEG(self):
        pass

    def plotEEGSpectrogram(self):
        pass

    def plotRaster(self):
        ### Function to create a raster plot after gathering the spikes ###
        print('Plotting raster ...')
        # Indices for each cell type based on their ID ranges
        pyso_indices = [i for i, val in enumerate(self.idVecAll) if val < self.Npyso]
        pydr_indices = [i for i, val in enumerate(self.idVecAll) if self.Npyso <= val < self.Npyso + self.Npydr]
        inh_indices = [i for i, val in enumerate(self.idVecAll) if self.Npyso + self.Npydr <= val < self.Npyso + self.Npydr + self.Ninh]
        if config.current_state_params['include_thalamus']:
            tc_indices = [i for i, val in enumerate(self.idVecAll) if self.Npyso + self.Npydr + self.Ninh <= val < self.Npyso + self.Npydr + self.Ninh + self.Ntc]
            re_indices = [i for i, val in enumerate(self.idVecAll) if self.Npyso + self.Npydr + self.Ninh + self.Ntc <= val < self.N]

        # Spike times and IDs for each cell type
        pyso_tVec = [self.tVecAll[i] for i in pyso_indices]
        pyso_id = [self.idVecAll[i] for i in pyso_indices]
        pydr_tVec = [self.tVecAll[i] for i in pydr_indices]
        pydr_id = [self.idVecAll[i] for i in pydr_indices]
        inh_tVec = [self.tVecAll[i] for i in inh_indices]
        inh_id = [self.idVecAll[i] for i in inh_indices]
        if config.current_state_params['include_thalamus']:
            tc_tVec = [self.tVecAll[i] for i in tc_indices]
            tc_id = [self.idVecAll[i] for i in tc_indices]
            re_tVec = [self.tVecAll[i] for i in re_indices]
            re_id = [self.idVecAll[i] for i in re_indices]

        import matplotlib.pyplot as plt

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
                (0, self.Npyso),  # PYso: 0 to Npyso-1
                (self.Npyso, self.Npyso + self.Npydr),  # PYdr: Npyso to Npyso+Npydr-1
                (self.Npyso + self.Npydr, self.Npyso + self.Npydr + self.Ninh),  # IN
                (self.Npyso + self.Npydr + self.Ninh, self.Npyso + self.Npydr + self.Ninh + self.Ntc),  # TC
                (self.Npyso + self.Npydr + self.Ninh + self.Ntc, self.N)  # RE
            ]
        else:
            id_ranges = [
                (0, self.Npyso),  # PYso: 0 to Npyso-1
                (self.Npyso, self.Npyso + self.Npydr),  # PYdr: Npyso to Npyso+Npydr-1
                (self.Npyso + self.Npydr, self.Npyso + self.Npydr + self.Ninh)  # IN
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

                # if label in ['PYso', 'PYdr']:
                #     # For PYso and PYdr, use scatter overlay for better visibility
                #     ax.scatter(tVec, ids, marker=',', s=0.5, color=color, alpha=0.3, label=label)
            else:
                # If no spikes, add a text annotation
                ax.text(0.5, 0.5, 'No spikes', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, color='gray')

            # Set y-axis limits and invert (0 at top, max at bottom)
            ax.set_ylim(id_min, id_max)
            ax.invert_yaxis()
            ax.set_ylabel(f'{label}')
            # ax.grid(True, linestyle='--', alpha=0.3)

        # Set x-axis label on the bottom subplot
        axes[-1].set_xlabel('Time (s)')

        # Set x-axis limits for all subplots
        if self.tVecAll:
            for ax in axes:
                ax.set_xlim(0, 1.01 * max(self.tVecAll) / 1000.0)

        # --- Save and close ---

        # Figure title - concise and informative
        # For publication: remove title (add as figure caption instead)
        # For presentations/internal: uncomment appropriate line
        # fig.suptitle('Propofol-Induced Network Activity', fontsize=10, y=0.995)
        
        # Save in multiple formats for journal submission
        if config.CONDITION == 'wake':
            base_filename = 'raster_plot_wake'
        elif config.CONDITION == 'propofol':
            base_filename = 'raster_plot_propofol'
        elif config.CONDITION == 'cortex_only_swo':
            base_filename = 'raster_plot_cortex_only_swo'
        elif config.CONDITION == 'disconnected_thalamus':
            base_filename = 'raster_plot_disconnected_thalamus'
        elif config.CONDITION == 'propofol_mostly_c':
            base_filename = 'raster_plot_propofol_mostly_c'
        elif config.CONDITION == 'propofol_mostly_i':
            base_filename = 'raster_plot_propofol_mostly_i'
        
        # High-resolution PNG for review/submission
        plt.savefig(f'{base_filename}.png', dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
        
        # # Vector format (EPS or PDF) - required by most journals
        # plt.savefig(f'{base_filename}.eps', bbox_inches='tight', 
        #         facecolor='white', edgecolor='none')
        
        # # TIFF format (some journals prefer this)
        # plt.savefig(f'{base_filename}.tiff', dpi=300, bbox_inches='tight',
        #         facecolor='white', edgecolor='none', pil_kwargs={"compression": "tiff_lzw"})
        
        plt.close()
        
        print(f'Publication-quality raster plot saved:')
        print(f'  - {base_filename}.png (300 DPI)')
        # print(f'  - {base_filename}.eps (vector)')
        # print(f'  - {base_filename}.tiff (300 DPI)')

    def printSpikeRates(self):
        ### Function to calculate the spike rate for each cell type ###
        print('Calculating and printing spike rates ...')
        
        # Indices for each cell type based on their ID ranges
        pyso_indices = [i for i, val in enumerate(self.idVecAll) if val < self.Npyso]
        pydr_indices = [i for i, val in enumerate(self.idVecAll) if self.Npyso <= val < self.Npyso + self.Npydr]
        inh_indices = [i for i, val in enumerate(self.idVecAll) if self.Npyso + self.Npydr <= val < self.Npyso + self.Npydr + self.Ninh]
        if config.current_state_params['include_thalamus']:
            tc_indices = [i for i, val in enumerate(self.idVecAll) if self.Npyso + self.Npydr + self.Ninh <= val < self.Npyso + self.Npydr + self.Ninh + self.Ntc]
            re_indices = [i for i, val in enumerate(self.idVecAll) if self.Npyso + self.Npydr + self.Ninh + self.Ntc <= val < self.N]

        # Spike times for each cell type
        pyso_tVec = [self.tVecAll[i] for i in pyso_indices]
        pydr_tVec = [self.tVecAll[i] for i in pydr_indices]
        inh_tVec = [self.tVecAll[i] for i in inh_indices]
        if config.current_state_params['include_thalamus']:
            tc_tVec = [self.tVecAll[i] for i in tc_indices]
            re_tVec = [self.tVecAll[i] for i in re_indices]

        # Number of cells in each type (from class attributes)
        if config.current_state_params['include_thalamus']:
            num_cells = {
                'PYso': self.Npyso,
                'PYdr': self.Npydr,
                'IN': self.Ninh,
                'TC': self.Ntc,
                'RE': self.N - (self.Npyso + self.Npydr + self.Ninh + self.Ntc)
            }
        else:
            num_cells = {
                'PYso': self.Npyso,
                'PYdr': self.Npydr,
                'IN': self.Ninh
            }

        # Determine simulation duration (in milliseconds)
        if self.tVecAll:
            sim_duration_ms = max(self.tVecAll)  # Maximum time in ms
        else:
            sim_duration_ms = 0  # Default to 0 if no spikes
        sim_duration_s = sim_duration_ms / 1000.0  # Convert to seconds

        # List of cell types and their spike times
        if config.current_state_params['include_thalamus']:
            cell_types = [
                ('PYso', pyso_tVec),
                ('PYdr', pydr_tVec),
                ('IN', inh_tVec),
                ('TC', tc_tVec),
                ('RE', re_tVec)
            ]
        else:
            cell_types = [
                ('PYso', pyso_tVec),
                ('PYdr', pydr_tVec),
                ('IN', inh_tVec)
            ]

        # Calculate and print spike rates
        for cell_type, tVec in cell_types:
            num_spikes = len(tVec)  # Total number of spikes for this cell type
            num_cells_type = num_cells[cell_type]
            
            if num_cells_type > 0 and sim_duration_s > 0:
                # Spike rate per cell in Hz (spikes/second)
                spike_rate = num_spikes / (num_cells_type * sim_duration_s)
            else:
                spike_rate = 0.0  # Avoid division by zero
            
            print(f"{cell_type} Spike Rate: {spike_rate:.2f} Hz ({num_spikes} spikes, {num_cells_type} cells, {sim_duration_ms:.2f} ms)")

    def saveData(self, output_base_dir=os.path.join('results', 'data')):
        if config.idhost == 0:
            if config.CONDITION == 'propofol_mostly_c' or config.CONDITION == 'propofol_mostly_i':
                print(f"Saving MINIMAL data for {config.CONDITION} | Stim: {config.stim_type} at {config.stim_start}ms")
            else:
                print(f'Saving data for condition: {config.CONDITION} ...')
            
            # 1. Base data package (Metadata and Spikes are always saved)
            dataSave = {
                'network_params': {
                    'Npyso': self.Npyso,
                    'Npydr': self.Npydr,
                    'Ninh': self.Ninh,
                    'Ntc': self.Ntc,
                    'Nre': self.Nre,
                    'condition': config.CONDITION,
                    'stim_type': config.stim_type,
                    'stim_start': config.stim_start,
                    'stim_dur': config.stim_dur,
                    'seed': config.randSeed, # Track the independent trial
                    'gAMPA_PY_PY': config.current_state_params['gAMPA_PYso_PYdr'],
                    'gAMPA_TC_PY': config.current_state_params['gAMPA_TC_PYdr'],
                    # Add flag status to metadata for tracking
                    'recorded_voltages': config.record_voltages 
                },
                'spikes': {
                    'tVec': self.tVecAll if hasattr(self, 'tVecAll') else [],
                    'idVec': self.idVecAll if hasattr(self, 'idVecAll') else []
                }
            }
            
            # 2. Conditionally compute EEG and append voltages
            if config.record_voltages:
                eeg_signal = np.zeros_like(self.t_vec_global)
                pydr_start = self.Npyso
                pydr_end = self.Npyso + self.Npydr
                
                for gid in range(pydr_start, pydr_end):
                    if gid in self.all_voltages:
                        eeg_signal += self.all_voltages[gid]
                
                dataSave['voltages'] = self.all_voltages
                dataSave['eeg'] = {
                    'time': self.t_vec_global,
                    'signal': eeg_signal
                }
                print(" -> Included full voltage arrays and computed EEG.")
            else:
                print(" -> Included minimal raster data only.")
            
            # 3. Create directory structure and save
            if config.CONDITION in ['propofol_mostly_c', 'propofol_mostly_i'] and config.stim_type == 'none':
                condition_dir = output_base_dir  # Save directly in base directory for minimal conditions
            else:
                condition_dir = os.path.join(output_base_dir, config.CONDITION)
            os.makedirs(condition_dir, exist_ok=True)
            
            # 4. Generate a unique filename using a timestamp
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            # Custom filename formatting based on the run type
            if config.args.table_name in ['table1', 'table2']:
                # Format: sim_table1_gTC_0.004_gPY_0.004_seed_1.pkl
                gTC_str = f"{config.current_state_params['gAMPA_TC_PYdr']*1e3:.4f}"
                gPY_str = f"{config.current_state_params['gAMPA_PYso_PYdr']*1e3:.4f}"
                filename = f"sim_{config.args.table_name}_gTC_{gTC_str}_gPY_{gPY_str}_seed_{config.randSeed}.pkl"
            elif config.CONDITION in ['propofol_mostly_c', 'propofol_mostly_i']:
                filename = f"sim_{config.CONDITION}_{config.stim_type}_{int(config.stim_start)}_{timestamp}.pkl"
            else:
                filename = f'sim_data_{config.CONDITION}_{timestamp}.pkl'

            filepath = os.path.join(condition_dir, filename)
            
            with open(filepath, 'wb') as f:
                pickle.dump(dataSave, f)
            print(f'Data successfully saved to {filepath}')