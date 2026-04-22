"""
Main simulation file for replicating 
"Rapid thalamocortical network switching mediated by cortical synchronization
underlies propofol-induced EEG signatures: a biophysical model," JNP, 2023.
"""
# --- import system libraries ---
import os
import sys

# --- import the required libraries ---
from neuron import h, coreneuron
from mpi4py import MPI
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import config

# ==========================================
# SMART CHECKPOINTING: Skip already completed runs
# ==========================================
if config.args.table_name in ['table1', 'table2']:
    # Reconstruct the exact filename this simulation would save
    gTC_str = f"{config.current_state_params['gAMPA_TC_PYdr']*1e3:.4f}"
    gPY_str = f"{config.current_state_params['gAMPA_PYso_PYdr']*1e3:.4f}"
    filename = f"sim_{config.args.table_name}_gTC_{gTC_str}_gPY_{gPY_str}_seed_{config.randSeed}.pkl"
    
    filepath = os.path.join('results', 'data', config.CONDITION, filename)
    
    # If the file is already on the hard drive, exit instantly AND cleanly
    if os.path.exists(filepath):
        if config.pc.id() == 0: # Only print from the master MPI node
            print(f"Skipping {filename}: Already exists!")
        
        config.pc.barrier() 
        config.pc.done()    
        h.quit()

h('load_file("stdgui.hoc")') # need this instead of import gui to get the simulation to be reproducible and not give an LFP flatline
from network import Net

# h('CVode[0].use_fast_imem(1)') #see use_fast_imem() at https://www.neuron.yale.edu/neuron/static/new_doc/simctrl/cvode.html

h.load_file("stdrun.hoc")

# ==========================================
# 1. ENABLE CORENEURON
# ==========================================
coreneuron.enable = config.args.coreneuron  # Set to True to enable CoreNEURON, False to run with standard NEURON (for testing/debugging/analysis)
coreneuron.gpu = config.args.gpu            # Set False since we built CoreNEURON for CPU
coreneuron.verbose = 0                      # Optional (To reduce the chatter)

# Optional: Print configuration to confirm
if config.idhost == 0:
    print(f"CoreNEURON Enabled: {coreneuron.enable}")
    print(f"GPU Mode: {coreneuron.gpu}")

def get_randomizer(id1, id2, seed):
    """
    Creates a Random123 object that is unique for:
    1. The specific cell/synapse (id1)
    2. The cell/synapse type (id2)
    3. The global simulation trial (seed)
    """
    r = h.Random()
    r.Random123(id1, id2, seed) 
    return r

def onerun(randSeed,Npyso,Npydr,Ninh,Ntc,Nre):
       
    # Random number generation and initialization   
    # h.Random().Random123_globalindex(randSeed) #this changes ALL Random123 streams
    # r = h.Random()
    
    # create network
    net = Net(Npyso,Npydr,Ninh,Ntc,Nre)

    # # Set up recording for cells you want to plot
    # if net.cells:  # Ensure cells exist
    #     net.cells[0].setRecording()  # Record traces for cell 0
         
    '''set up custom initialization to set initial values of all cells'''
    def set_initial_conditions(randSeed):
        if config.idhost == 0: print(f"Initializing network with Global Seed: {randSeed}")

        #loop through all PYso cells and set their initial conditions uniquely
        for pyso_gid in net.pyso_gidList:
            if config.pc.gid_exists(pyso_gid): # CRITICAL CHECK
                # r.Random123(pyso_gid,0,0) #set stream of random numbers; first argument is gid, make second argument 0
                # Use helper: ID1=gid, ID2=0 (type), ID3=randSeed
                r = get_randomizer(pyso_gid, 0, randSeed)
                config.pc.gid2cell(pyso_gid).soma(0.5).hNa_PYso = 0.7 + (r.uniform(0,1)) * 0.1
                config.pc.gid2cell(pyso_gid).soma(0.5).nK_PYso = 0.05 + (r.uniform(0,1)) * 0.05
                config.pc.gid2cell(pyso_gid).soma(0.5).hA_PYso = 0.1 + (r.uniform(0,1)) * 0.1
                config.pc.gid2cell(pyso_gid).soma(0.5).mKS_PYso = 0.005 + (r.uniform(0,1)) * 0.001
                config.pc.gid2cell(pyso_gid).soma.v = -68 + (r.uniform(0,1)) * 20

        #loop through all PYdr cells and set their initial conditions uniquely
        for pydr_gid in net.pydr_gidList:
            if config.pc.gid_exists(pydr_gid): # CRITICAL CHECK
                # r.Random123(pydr_gid,1,0) #set stream of random numbers; first argument is gid, make second argument 1 
                r = get_randomizer(pydr_gid, 1, randSeed)
                config.pc.gid2cell(pydr_gid).soma(0.5).CaBuffer_PYdr = 0.001 + (r.uniform(0,1)) * 0.01
                config.pc.gid2cell(pydr_gid).soma.v = -68 + (r.uniform(0,1)) * 20

        #loop through all INH cells and set their initial conditions uniquely
        for inh_gid in net.inh_gidList:
            if config.pc.gid_exists(inh_gid): # CRITICAL CHECK
                # r.Random123(inh_gid,2,0) #set stream of random numbers; first argument is gid, make second argument 2
                r = get_randomizer(inh_gid, 2, randSeed)
                config.pc.gid2cell(inh_gid).soma(0.5).hNa_IN = 0.7 + (r.uniform(0,1)) * 0.1
                config.pc.gid2cell(inh_gid).soma(0.5).nK_IN = 0.13 + (r.uniform(0,1)) * 0.01
                config.pc.gid2cell(inh_gid).soma.v = -68 + (r.uniform(0,1)) * 20

        if config.current_state_params['include_thalamus']:
            #loop through all TC cells and set their initial conditions uniquely
            for tc_gid in net.tc_gidList:
                if config.pc.gid_exists(tc_gid): # CRITICAL CHECK
                    # r.Random123(tc_gid,3,0) #set stream of random numbers; first argument is gid, make second argument 3
                    r = get_randomizer(tc_gid, 3, randSeed)
                    config.pc.gid2cell(tc_gid).soma(0.5).cai_TC = 0.0003 + (r.uniform(0,1)) * 0.00001
                    config.pc.gid2cell(tc_gid).soma(0.5).mNa_TC = 0.00007 + (r.uniform(0,1)) * 0.00001
                    config.pc.gid2cell(tc_gid).soma(0.5).hNa_TC = 0.8 + (r.uniform(0,1)) * 0.1
                    config.pc.gid2cell(tc_gid).soma(0.5).nK_TC = 0.000250 + (r.uniform(0,1)) * 0.00001
                    config.pc.gid2cell(tc_gid).soma(0.5).hT_TC = 0.01 + (r.uniform(0,1)) * 0.005
                    config.pc.gid2cell(tc_gid).soma(0.5).Open_TC = 0.05 + (r.uniform(0,1)) * 0.01
                    config.pc.gid2cell(tc_gid).soma(0.5).Pone_TC = 0.06 + (r.uniform(0,1)) * 0.01
                    config.pc.gid2cell(tc_gid).soma(0.5).OpenLocked_TC = 0.55 + (r.uniform(0,1)) * 0.01
                    config.pc.gid2cell(tc_gid).soma.v = -68 + (r.uniform(0,1)) * 20           

            #loop through all RE cells and set their initial conditions uniquely
            for re_gid in net.re_gidList:
                if config.pc.gid_exists(re_gid): # CRITICAL CHECK
                    # r.Random123(re_gid,4,0) #set stream of random numbers; first argument is gid, make second argument 4
                    r = get_randomizer(re_gid, 4, randSeed)
                    config.pc.gid2cell(re_gid).soma(0.5).mNa_RE = 0.00002 + (r.uniform(0,1)) * 0.00001
                    config.pc.gid2cell(re_gid).soma(0.5).hNa_RE = 0.8 + (r.uniform(0,1)) * 0.1
                    config.pc.gid2cell(re_gid).soma(0.5).nK_RE = 0.00015 + (r.uniform(0,1)) * 0.00001
                    config.pc.gid2cell(re_gid).soma(0.5).mT_RE = 0.01 + (r.uniform(0,1)) * 0.01
                    config.pc.gid2cell(re_gid).soma(0.5).hT_RE = 0.6 + (r.uniform(0,1)) * 0.01
                    config.pc.gid2cell(re_gid).soma.v = -68 + (r.uniform(0,1)) * 20

        # --- Direct compartmental connections ---

        for i, syn in enumerate(net.PYso_PYdr_ICOM_synapses):
            # r.Random123(i, 5, 0)
            r = get_randomizer(i, 5, randSeed)
            # syn.gCOM = 0.005 + (r.uniform(-1,1)) * 0.000286 # S/cm2
            syn.gCOM = r.normal(0.005, 0.000286**2)  # S/cm2
            # syn.gCOM = 0.005 # S/cm2 (no variability)

        for i, syn in enumerate(net.PYdr_PYso_ICOM_synapses):
            # r.Random123(i, 5, 0)
            r = get_randomizer(i, 5, randSeed)
            # syn.gCOM = 0.011667 + (r.uniform(-1,1)) * 0.000667 # S/cm2
            syn.gCOM = r.normal(0.011667, 0.000667**2)  # S/cm2
            # syn.gCOM = 0.011667 # S/cm2 (no variability)

        for i, syn in enumerate(net.PYdr_PYso_IKNa_synapses):
            # r.Random123(i, 6, 0)
            r = get_randomizer(i, 6, randSeed)
            syn.concNa = 12 + (r.uniform(0,1)) * 3
            syn.hNalocal = 0.5 + (r.uniform(0,1)) * 0.1

        # --- Intracortical synapses ---

        for i, syn in enumerate(net.AMPA_PYso_PYdr_synapses):
            # r.Random123(i, 7, 0)
            r = get_randomizer(i, 7, randSeed)
            syn.sAMPA = 0.0 + 0.1 * r.uniform(0, 1)
            syn.res_AMPA = 0.8 - 0.1 * r.uniform(0, 1)

        for i, syn in enumerate(net.NMDA_PYso_PYdr_synapses):
            # r.Random123(i, 8, 0)
            r = get_randomizer(i, 8, randSeed)
            syn.sNMDA = 0.0 + 0.1 * r.uniform(0, 1)
            syn.xNMDA = 0.0 + 0.1 * r.uniform(0, 1)
            syn.res_NMDA = 0.8 - 0.1 * r.uniform(0, 1)

        for i, syn in enumerate(net.AMPA_PYso_IN_synapses):
            # r.Random123(i, 9, 0)
            r = get_randomizer(i, 9, randSeed)
            syn.sAMPA = 0.0 + 0.1 * r.uniform(0, 1)
            syn.res_AMPA = 0.8 - 0.1 * r.uniform(0, 1)

        for i, syn in enumerate(net.NMDA_PYso_IN_synapses):
            # r.Random123(i, 10, 0)
            r = get_randomizer(i, 10, randSeed)
            syn.sNMDA = 0.0 + 0.1 * r.uniform(0, 1)
            syn.xNMDA = 0.0 + 0.1 * r.uniform(0, 1)
            syn.res_NMDA = 0.8 - 0.1 * r.uniform(0, 1)

        for i, syn in enumerate(net.GABAA_IN_PYso_synapses):
            # r.Random123(i, 11, 0)
            r = get_randomizer(i, 11, randSeed)
            syn.sGABAA = 0.0 + 0.1 * r.uniform(0, 1)
            syn.res_GABAA = 0.8 - 0.1 * r.uniform(0, 1)

        for i, syn in enumerate(net.GABAA_IN_IN_synapses):
            # r.Random123(i, 12, 0)
            r = get_randomizer(i, 12, randSeed)
            syn.sGABAA = 0.0 + 0.1 * r.uniform(0, 1)
            syn.res_GABAA = 0.8 - 0.1 * r.uniform(0, 1)

        # --- Thalamic synapses ---

        if config.current_state_params['include_thalamus']:
            for i, syn in enumerate(net.AMPA_TC_RE_synapses):
                # r.Random123(i, 13, 0)
                r = get_randomizer(i, 13, randSeed)
                syn.sAMPA = 0.0 + 0.1 * r.uniform(0, 1)

            for i, syn in enumerate(net.GABAA_RE_RE_synapses):
                # r.Random123(i, 14, 0)
                r = get_randomizer(i, 14, randSeed)
                syn.sGABAA = 0.0 + 0.1 * r.uniform(0, 1)

            for i, syn in enumerate(net.GABAA_RE_TC_synapses):
                # r.Random123(i, 15, 0)
                r = get_randomizer(i, 15, randSeed)
                syn.sGABAA = 0.0 + 0.1 * r.uniform(0, 1)

            for i, syn in enumerate(net.GABAB_RE_TC_synapses):
                # r.Random123(i, 16, 0)
                r = get_randomizer(i, 16, randSeed)
                syn.rGABAB = 0.0 + 0.1 * r.uniform(0, 1)
                syn.sGABAB = 0.0 + 0.1 * r.uniform(0, 1)
        
        # --- Thalamocortical synapses ---

        if config.current_state_params['include_thalamocortical_connections']:
            for i, syn in enumerate(net.AMPA_TC_PYdr_synapses):
                # r.Random123(i, 17, 0)
                r = get_randomizer(i, 17, randSeed)
                syn.sAMPA = 0.0 + 0.1 * r.uniform(0, 1)

            for i, syn in enumerate(net.AMPA_TC_IN_synapses):
                # r.Random123(i, 18, 0)
                r = get_randomizer(i, 18, randSeed)
                syn.sAMPA = 0.0 + 0.1 * r.uniform(0, 1)

            for i, syn in enumerate(net.AMPA_PYso_TC_synapses):
                # r.Random123(i, 19, 0)
                r = get_randomizer(i, 19, randSeed)
                syn.sAMPA = 0.0 + 0.1 * r.uniform(0, 1)

            for i, syn in enumerate(net.AMPA_PYso_RE_synapses):
                # r.Random123(i, 20, 0)
                r = get_randomizer(i, 20, randSeed)
                syn.sAMPA = 0.0 + 0.1 * r.uniform(0, 1)
        
        if config.idhost == 0:
            print("Custom initialization complete.") 
          
    # Wrap the initialization in a lambda or partial to pass the seed
    # The FInitializeHandler needs a function with NO arguments, 
    # so we create a simple wrapper function here.
    def init_wrapper():
        set_initial_conditions(randSeed)

    # run sim and gather spikes
    h.dt = config.time_step #set the time step as set in the config file
    h.steps_per_ms = 1/h.dt #set the number of steps per ms
    config.pc.set_maxstep(10)  # ms, maximum step size for parallel simulation; see https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/network/parcon.html#ParallelContext.set_maxstep
    # see https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/network/parcon.html#ParallelContext.set_maxstep
    # as well as section 2.4 of "Simulation Neurotechnologies for Advancing Brain Research: 
    # Parallelizing Large Networks in NEURON" (Lytton et. al, 2016)
    
    # Register the wrapper
    fih = h.FInitializeHandler(init_wrapper)

    h.finitialize(-68) #set initial voltages of all cells to -68 mV or default values
    # h.stdinit()
    # h.init()  # Runs INITIAL blocks with current v values
    
    # ==========================================
    # 2. RUN SIMULATION (Optimized with Progress Bar)
    # ==========================================
    if config.idhost == 0: 
        print(f'Running sim with CoreNEURON: {coreneuron.enable}...')
        startTime = datetime.now() 

    run_duration = config.duration
    report_interval = 1000.0  # Report progress every 1000 ms simulation time
    current_time = 0.0
    
    # Loop over large chunks of time to allow CoreNEURON to optimize,
    # but interrupt occasionally to print progress.
    while current_time < run_duration:
        next_stop = min(current_time + report_interval, run_duration)
        
        # Run CoreNEURON for this chunk
        config.pc.psolve(next_stop)
        
        current_time = next_stop
        
        # Print Progress (Rank 0 only)
        if config.idhost == 0:
            percent = (current_time / run_duration) * 100
            # \r overwrites the line, keeping the terminal clean
            print(f"  [Progress] {current_time:.0f} ms / {run_duration:.0f} ms ({percent:.1f}%)")
    
    # ==========================================
    # 3. GATHER RESULTS
    # ==========================================
    
    if config.idhost == 0:
        runTime = (datetime.now() - startTime).total_seconds()
        print(f"Run time for {int(config.duration/1000.0)} sec sim = {runTime:.2f} sec")
    
    net.gatherSpikes()

    # Only gather distributed voltages if the flag was passed
    if config.record_voltages:
        net.gatherVoltages()
    
    if config.idhost == 0: 
        # Save raster to file
        with open(f"raster_nhost={config.nhost}.txt", 'w') as f:
            for i in range(len(net.tVecAll)):
                f.write(f"{net.tVecAll[i]:.3f}  {net.idVecAll[i]:.0f}\n")
        
        # Plotting
        if config.CONDITION != 'tables_propofol':
            if len(net.tVecAll) > 0:
                net.plotRaster()
                net.printSpikeRates()
            else:
                print("No cells spiked.")
        
        net.saveData()

        # # --- Test Plotting of Voltage Traces ---
        # if net.cells:  # Check if cells exist
        #     net.cells[0].plotTraces()  # Plot voltage traces for cell 0
        
    del net  # free memory

# Run the simulation
onerun(config.randSeed,config.Npyso,config.Npydr,config.Ninh,config.Ntc,config.Nre)

# ==========================================
# GRACEFUL MPI SHUTDOWN
# ==========================================
config.pc.barrier() # Wait for all MPI nodes to catch up (Rank 5 waits for Rank 0 to finish saving)
config.pc.done()    # Tells NEURON to cleanly wrap up MPI
h.quit()            # Cleanly exit the NEURON environment
