#!/bin/bash

# Trap Ctrl+C (SIGINT) and exit the entire script immediately
trap "echo 'Batch execution cancelled by user.'; exit 1" SIGINT

# Define simulation boundaries (in milliseconds)
START_TIME=7500
END_TIME=14900
STEP_SIZE=100

echo "Starting exhaustive batch simulations for Figures 4 and 5..."

# Loop from 0 to 14900 in steps of 100
for stim_start in $(seq $START_TIME $STEP_SIZE $END_TIME); do
    echo "========================================================="
    echo "Running stimulus sweep for interval starting at: ${stim_start} ms"
    echo "========================================================="
    
    # Figure 4 Experiments: C-State to I-State (Baseline is mostly C-state)
    mpiexec -n 8 nrniv -mpi -python SoplataNet.py --coreneuron --condition propofol_mostly_c --stim_type sync --stim_start $stim_start
    mpiexec -n 8 nrniv -mpi -python SoplataNet.py --coreneuron --condition propofol_mostly_c --stim_type desync --stim_start $stim_start
    
    # Figure 5 Experiments: I-State to C-State (Baseline is mostly I-state)
    mpiexec -n 8 nrniv -mpi -python SoplataNet.py --coreneuron --condition propofol_mostly_i --stim_type sync --stim_start $stim_start
    mpiexec -n 8 nrniv -mpi -python SoplataNet.py --coreneuron --condition propofol_mostly_i --stim_type desync --stim_start $stim_start

done

echo "All 300 batch simulations completed successfully!"
