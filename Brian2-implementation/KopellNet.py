"""
KopellNet.py — Main Entry Point
=====================================
Builds the thalamocortical network, runs the simulation, and then calls
the standard analysis pipeline (raster plot, spike rates, data save).

This script mirrors the role of SoplataNet.py in the NEURON implementation.
All configuration is driven by config.py and its command-line arguments, so
this file should not normally need to be edited.

Typical usage
-------------
# Default propofol condition (30 s, 100 PYso/PYdr, 20 IN/TC/RE):
python KopellNet.py --condition propofol --seed 101

# Wake-state with voltage traces recorded:
python KopellNet.py --condition wake --record_voltages

# Tables 1 & 2 batch sweep (override conductances):
python KopellNet.py --condition tables_propofol \\
    --table_name table1 --gAMPA_PY_PY 0.005e-3 --gAMPA_TC_PY 0.005e-3 --seed 42

# Cortex-only SWO, 2 OpenMP threads:
python KopellNet.py --condition cortex_only_swo --openmp_threads 2
"""

import time
import numpy as np

# config.py must be imported first — it sets up the Brian2 device & clock
import config
from network import Net


def main():
    # -------------------------------------------------------------------------
    # Reproducibility
    # -------------------------------------------------------------------------
    np.random.seed(config.randSeed)

    # -------------------------------------------------------------------------
    # Build and run
    # -------------------------------------------------------------------------
    t_start = time.time()

    net = Net()
    net.run()

    t_elapsed = time.time() - t_start
    print(f'Total wall-clock time: {t_elapsed:.1f} s  '
          f'({t_elapsed/60:.2f} min)')

    # -------------------------------------------------------------------------
    # Analysis pipeline  (mirrors NEURON SoplataNet.py post-run calls)
    # -------------------------------------------------------------------------
    net.printSpikeRates()
    net.plotRaster(output_dir='results')

    if config.record_voltages:
        # Plot a sample voltage trace for the first PYso and first TC cell
        net.plotVoltageTrace(cell_type='PYso', cell_idx=0, output_dir='results')
        if config.current_state_params['include_thalamus']:
            net.plotVoltageTrace(cell_type='TC', cell_idx=0, output_dir='results')

    net.saveData(output_base_dir='results/data')


if __name__ == '__main__':
    main()
