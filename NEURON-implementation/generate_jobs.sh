#!/bin/bash

# Clear any existing file
> jobs_list.txt

# Use your explicit conda python path to prevent environment errors on the remote machines
PYTHON_CMD="mpiexec -n 4 --bind-to none nrniv -mpi -python SoplataNet.py"

echo "Generating Table 1 commands..."
for g in 0.002e-3 0.003e-3 0.004e-3 0.005e-3 0.006e-3 0.007e-3 0.008e-3 0.009e-3 0.01e-3 0.012e-3; do
    for seed in $(seq 1 70); do
        echo "$PYTHON_CMD --coreneuron --condition tables_propofol --table_name table1 --gAMPA_PY_PY $g --gAMPA_TC_PY $g --seed $seed" >> jobs_list.txt
    done
done

echo "Generating Table 2 commands..."
for gTC in 0.004e-3 0.006e-3 0.008e-3 0.01e-3; do
    for gPY in 0.004e-3 0.006e-3 0.008e-3 0.01e-3; do
        for seed in $(seq 1 40); do
            echo "$PYTHON_CMD --coreneuron --condition tables_propofol --table_name table2 --gAMPA_TC_PY $gTC --gAMPA_PY_PY $gPY --seed $seed" >> jobs_list.txt
        done
    done
done

echo "Successfully generated 1,340 commands in jobs_list.txt"
