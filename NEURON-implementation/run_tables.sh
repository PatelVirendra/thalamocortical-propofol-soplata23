#!/bin/bash

# Define the number of concurrent simulations.
# Formula: (Total Physical Cores / 4) - 1 buffer space
# Assuming that the physical cores of your system are in multiple of 4 
CONCURRENT_JOBS=11 # For a 48 core cloud VM

echo "Starting 1,340 simulations via GNU Parallel using $CONCURRENT_JOBS concurrent workers..."

# Read the jobs_list and feed it to the parallel dispatcher
# The --bar flag creates the live progress UI
parallel -j $CONCURRENT_JOBS --bar --joblog parallel_task_log.txt < jobs_list.txt

echo "All batch runs completed on the Cloud VM!"
