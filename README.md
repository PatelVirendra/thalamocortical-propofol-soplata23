# Biophysical Modeling of Thalamocortical Dynamics Under Propofol Anesthesia

This repository contains large-scale biophysical thalamocortical network models implemented in both **NEURON (with CoreNEURON)** and **Brian2**. 

The models are designed to replicate the oscillatory dynamics and network state transitions observed under propofol-induced unconsciousness, as originally described in:

> **Rapid thalamocortical network switching mediated by cortical synchronization underlies propofol-induced EEG signatures: a biophysical model.** > *Austin E. Soplata, Elie Adam, Emery N. Brown, Patrick L. Purdon, Michelle M. McCarthy, and Nancy Kopell. (Journal of Neurophysiology, 2023).*

## 🧠 Scientific Background

Propofol-mediated unconsciousness is characterized by strong alpha/low-beta (8–20 Hz) and slow (0.1–1 Hz) oscillations in the EEG. As the anesthetic dose increases, the spectral signature changes, indicating deeper levels of unconsciousness. 

This biophysical model explains these transitions by showing that the thalamocortical network fluctuates between two mutually exclusive states:
1.  **The C-state (Continuous):** Characterized by continuous alpha/low-beta spiking in the thalamus, low cortical synchrony, and positive corticothalamic feedback.
2.  **The I-state (Interrupted):** Characterized by highly synchronized cortical slow waves that periodically interrupt thalamic spiking, forcing a homeostatic feedback loop. 

As propofol dose increases (modeled via changes to thalamocortical and intracortical synaptic conductances driven by lowered cholinergic tone), the network spends increasingly more time in the I-state, perfectly recapitulating the clinical EEG phenomenology.

---

## 📂 Repository Structure

```text
├── NEURON_Engine/          # High-performance MPI/CoreNEURON implementation
│   ├── mod-files/          # NMODL mechanism files (.mod)
│   ├── config.py           # Global parameters and state switches
│   ├── cells.py            # Cell class definitions (PYso, PYdr, IN, TC, RE)
│   ├── netcon.py           # Nearest Neighbors Connectivity function
│   ├── network.py          # Synapses, Spike and EEG gathering, and Data saving
│   ├── SoplataNet.py       # Main execution script for the NEURON model
│   ├── generate_jobs.sh    # Bash script for generating simulation jobs for Table 1 & 2 from the paper
│   ├── run_experiments.sh  # Bash script for running simulations for generating Fig. 4 & 5 from the paper
│   └── run_tables.sh       # Bash script for running simulations jobs for Table 1 & 2 from the paper using GNU Parallel
│
├── Brian2_Engine/          # Brian2 implementation for rapid prototyping
│   ├── config.py           # Global parameters and state switches
│   ├── utils.py            # Randomization function for reproducibility (mimics Random123 from NEURON)
│   ├── cells.py            # Cell class definitions (PYso, PYdr, IN, TC, RE)
│   ├── netcon.py           # Nearest Neighbors Connectivity function
│   ├── synapses.py         # Synaptic connectivity logic
│   ├── network.py          # Builds and assemble the model. Includes spike and EEG gathering, and data saving functions
│   └── KopellNet.py        # Main execution script for the Brian2 model
│
├── Analysis_Pipeline/      # Shared scripts for generating paper figures
│   ├── analysis.py               # Time-series, spectrograms, and raster plotting
│   ├── analyze_perturbations.py  # HMM-based analysis of Sync/Desync stimuli (Fig 4 & 5)
│   └── analyze_tables.py         # HMM-based parameter sweep aggregations (Tables 1 & 2)
│
├── environment.yml         # Conda environment setup file
└── README.md
