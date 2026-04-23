# Biophysical Modeling of Thalamocortical Dynamics Under Propofol Anesthesia

This project develops a large-scale biophysical thalamocortical network model replicating the oscillatory dynamics and network state transitions observed under propofol-induced unconsciousness, as described in:

> **Rapid thalamocortical network switching mediated by cortical synchronization underlies propofol-induced EEG signatures: a biophysical model.**
> *Austin E. Soplata, Elie Adam, Emery N. Brown, Patrick L. Purdon, Michelle M. McCarthy, and Nancy Kopell. (Journal of Neurophysiology, 2023).*

Two implementations are provided: a high-performance **NEURON + CoreNEURON** engine for HPC environments, and a readable **Brian2** engine for rapid prototyping.

---

## Table of Contents

- [Scientific Background](#scientific-background)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
  - [Building NEURON with CoreNEURON from Source](#building-neuron-with-coreneuron-from-source)
  - [Compiling NMODL Mechanism Files](#compiling-nmodl-mechanism-files)
- [Running the Simulations](#running-the-simulations)
- [Generating Figures & Analysis](#generating-figures--analysis)
- [HPC Optimizations (CoreNEURON)](#hpc-optimizations-coreneuron)
- [Key Observations](#key-observations)
- [Key Issues in the Original DynaSim Code](#key-issues-in-the-original-dynasim-code)
- [Credits](#credits)
- [Author](#author)
- [License](#license)

---

## Scientific Background

Propofol-mediated unconsciousness is characterized by strong alpha/low-beta (8–20 Hz) and slow (0.1–1 Hz) oscillations in the EEG. As the anesthetic dose increases, the spectral signature shifts, indicating deeper levels of unconsciousness.

This model explains these transitions by showing that the thalamocortical network fluctuates between two mutually exclusive states:

1. **C-state (Continuous):** Characterized by continuous alpha/low-beta spiking in the thalamus, low cortical synchrony, and positive corticothalamic feedback.
2. **I-state (Interrupted):** Characterized by highly synchronized cortical slow waves that periodically interrupt thalamic spiking, forcing a homeostatic feedback loop.

As propofol dose increases — modeled by changes in thalamocortical and intracortical synaptic conductances driven by reduced cholinergic tone — the network spends more time in the I-state, recapitulating the clinical EEG phenomenology.

---

## Repository Structure

```text
├── NEURON-implementation/            # High-performance MPI/CoreNEURON implementation
│   ├── mod-files/                    # NMODL mechanism files (.mod)
│   ├── config.py                     # Global parameters and state switches
│   ├── cells.py                      # Cell class definitions (PYso, PYdr, IN, TC, RE)
│   ├── netcon.py                     # Nearest-neighbors connectivity function
│   ├── network.py                    # Synapses, spike/EEG gathering, and data saving
│   ├── SoplataNet.py                 # Main execution script
│   ├── generate_jobs.sh              # Generates simulation jobs for Tables 1 & 2
│   ├── run_experiments.sh            # Runs simulations for Figures 4 & 5
│   └── run_tables.sh                 # Runs Table 1 & 2 jobs via GNU Parallel
│
├── Brian2-implementation/            # Brian2 implementation for rapid prototyping
│   ├── config.py                     # Global parameters and state switches
│   ├── utils.py                      # Randomization function (mimics NEURON's Random123)
│   ├── cells.py                      # Cell class definitions (PYso, PYdr, IN, TC, RE)
│   ├── netcon.py                     # Nearest-neighbors connectivity function
│   ├── synapses.py                   # Synaptic connectivity logic
│   ├── network.py                    # Builds and assembles the model; data saving
│   └── KopellNet.py                  # Main execution script
│
├── Analysis_Pipeline/                # Shared scripts for generating paper figures
│   ├── analysis.py                   # Time-series, spectrograms, and raster plots
│   ├── analyze_perturbations.py      # HMM analysis of Sync/Desync stimuli (Figs. 4 & 5)
│   └── analyze_tables.py             # HMM-based parameter sweep aggregation (Tables 1 & 2)
│
├── environment.yml                   # Conda environment specification
└── README.md
```

| Feature | NEURON + CoreNEURON | Brian2 |
|---|---|---|
| **Best for** | Full-scale HPC runs | Rapid prototyping |
| **Parallelization** | MPI (`ParallelContext`) + OpenMP | OpenMP only |
| **Backend** | CoreNEURON (vectorized) | Pure Python/Brian2 |
| **Determinism** | Bitwise identical across same `nhost` | Standard |

---

## Installation

### 1. Set Up the Conda Environment

An `environment.yml` is provided to replicate the exact Python environment, including dependencies for simulation, signal processing, and HMM-based state classification.

```bash
conda env create -f environment.yml
conda activate anesthesia
```

> Brian2 requires no separate installation step — it is included in `environment.yml`. To install it standalone: `pip install brian2`.

---

### Building NEURON with CoreNEURON from Source

The simulation requires a **custom NEURON build** with CoreNEURON, MPI, and math optimizations enabled. The standard `pip install neuron` distribution does not include CoreNEURON support. Follow the steps below to compile from source.

#### Prerequisites

Ensure the following system packages are available before starting:

```bash
# On Debian/Ubuntu systems
sudo apt install -y git cmake build-essential gcc g++ \
    libncurses-dev libopenmpi-dev openmpi-bin \
    python3-dev flex bison wget
```

Also, install the required Python build tools into your conda environment:

```bash
conda activate anesthesia
pip install jinja2 setuptools wheel pyyaml cython
```

#### Step 1: Clone the NEURON Source

```bash
cd $HOME
git clone https://github.com/neuronsimulator/nrn.git
cd nrn
mkdir build
cd build
```

#### Step 2: Configure with CMake

Run CMake with the following flags. This enables CoreNEURON, MPI parallelization, SIMD math optimizations, and dynamic Python linking, while disabling the GUI (Interviews) and 3D morphology (RX3D) components that are not needed for this model.

```bash
cmake .. \
    -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
    -DNRN_ENABLE_CORENEURON=ON \
    -DNRN_ENABLE_MPI=ON \
    -DNRN_ENABLE_MATH_OPT=ON \
    -DNRN_ENABLE_PYTHON=ON \
    -DNRN_ENABLE_PYTHON_DYNAMIC=ON \
    -DPYTHON_EXECUTABLE=$(which python) \
    -DNRN_ENABLE_INTERVIEWS=OFF \
    -DNRN_ENABLE_RX3D=OFF
```

> **Important:** Always activate the conda environment (`conda activate anesthesia`) before running CMake so that `$(which python)` resolves to the correct interpreter. If CMake fails partway through, clear the build directory with `rm -rf *` before re-running.

#### Step 3: Compile and Install

```bash
# Adjust -j to your CPU core count (e.g., -j 8 for 8 cores)
make -j 8

# Install into the active conda environment ($CONDA_PREFIX)
make install
```

#### Step 4: Configure Environment Variables

NEURON's NMODL compiler requires a runtime variable pointing to the Python shared library. Add the following to your `~/.bashrc` (or `~/.zshrc`) so it persists across sessions:

```bash
# Ensure nrnivmodl can locate the Python shared library at runtime
echo 'export NMODL_PYLIB=$(python -c "import sysconfig, os; print(os.path.join(sysconfig.get_config_var('"'"'LIBDIR'"'"'), sysconfig.get_config_var('"'"'INSTSONAME'"'"')))")' >> ~/.bashrc

# Point NMODL to the conda prefix
echo 'export NMODLHOME=$CONDA_PREFIX' >> ~/.bashrc

source ~/.bashrc
conda activate anesthesia
```

#### Step 5: Verify the Installation

```bash
# Should print "MPI Threads: 1" when run outside of mpiexec
python -c "from neuron import h; pc = h.ParallelContext(); print(f'MPI Threads: {pc.nhost()}')"

# Full MPI test across 4 ranks
mpiexec -n 4 python -c "from neuron import h; pc = h.ParallelContext(); print(f'Rank {pc.id()} of {pc.nhost()}')"
```

> If using VSCode, update the Python interpreter path in the `.env` file located in the `NEURON-implementation/` directory to point to your conda environment.

---

### Compiling NMODL Mechanism Files

With NEURON installed, compile the `.mod` files in the `NEURON-implementation/` directory. This generates the `x86_64/` shared library that NEURON loads at runtime.

```bash
cd NEURON-implementation

# If recompiling after code changes, remove old build artifacts first
rm -rf x86_64

# Compile for CoreNEURON (generates CoreNEURON-optimized binaries)
nrnivmodl -coreneuron mod-files
```

A successful compilation produces an `x86_64/` directory containing `special` and `libnrnmech.so`.

> To use standard NEURON without CoreNEURON, compile with `nrnivmodl mod-files` (no `-coreneuron` flag) and omit `--coreneuron` from the run command below.

---

## Running the Simulations

### NEURON + CoreNEURON (Recommended for full-scale runs)

Simulation conditions (Wake, Propofol, Cortex-Only-SWO, etc.), key parameters, and stimulation type are all controlled in `config.py`. You can just browse through it for the full list of available options. Simulation duration is also set there.

```bash
cd NEURON-implementation

# Adjust -n based on your machine's core count
mpiexec -n 8 nrniv -mpi -python SoplataNet.py --coreneuron
```

> **Note:** The `--coreneuron` flag and voltage recording (`--record-voltages`) are mutually exclusive. CoreNEURON is optimized for memory efficiency and only transfers spike data back to NEURON. To record voltages, omit `--coreneuron` and recompile with `nrnivmodl mod-files`.

### Brian2 (Recommended for prototyping)

```bash
cd Brian2-implementation
python KopellNet.py
```

The default condition is propofol. See `config.py` and `KopellNet.py` for the full list of conditions and configuration options.

---

## Generating Figures & Analysis

Both engines output spike times, voltages, EEG, and micro-variables consumed by the shared scripts in `Analysis_Pipeline/`, replicating the exact figures and tables from Soplata et al. (2023).

### Step 1: Collect Simulation Data

Run the simulation with the `--condition` flag to produce data for conditions other than the default (propofol), such as Wake or Cortex-Only-SWO. Output `.pkl` files are saved to `results/data/<CONDITION>/`. A raster plot for each condition is also saved to the root directory.

### Step 2: Run Analysis Scripts

All analysis scripts default to reading data from `NEURON-implementation/`. To analyze Brian2 output, append `--engine Brian2-implementation` to any command below.

#### `analysis.py` → Figures 2 & 3

Applies Butterworth bandpass filtering, computes Welch's Power Spectral Density, and generates multitaper spectrograms. Output figures are saved to `results/figures/<CONDITION>/`.

- **Outputs:** Full EEG overview, zoomed raster plot, spectrogram (Fig. 2), and I-state/C-state micro-variable traces, raster plots, and filtered EEG (Fig. 3).

```bash
python Analysis_Pipeline/analysis.py
```

#### `analyze_perturbations.py` → Figures 4 & 5

Replicates stimulus perturbation experiments. Computes the delta/beta ratio from pyramidal spikes and uses a pre-trained Gaussian HMM to classify each time window as C-state or I-state.

- **Outputs:** Stacked bar charts showing the proportion of post-stimulus time converted from C-to-I or I-to-C after synchronizing/desynchronizing inputs.

> To generate perturbation data, use the `--stim_type` and `--stim_start` arguments in `SoplataNet.py` to inject a 100 ms `IClamp` during an ongoing C-state or I-state. The helper script `run_experiments.sh` automates this across conditions.

```bash
python Analysis_Pipeline/analyze_perturbations.py
```

#### `analyze_tables.py` → Tables 1 & 2

Scales HMM analysis across large parameter sweeps, scoring hundreds of simulation files over a 2D grid of intracortical (`gPY`) and thalamocortical (`gTC`) conductances.

- **Outputs:** Formatted terminal tables showing mean ± std of the percentage of time the network spends in C-state vs. I-state across random seeds.

> Run `generate_jobs.sh` first, then use `run_tables.sh` with GNU Parallel to batch-generate the required simulation files.

```bash
python Analysis_Pipeline/analyze_tables.py
```

> **Note on automation:** The original DynaSim code generates simulation files but relies on manual scoring of I-state and C-state percentages. Here, both states are scored automatically using an HMM trained and parameter-tuned on a standard propofol simulation from the paper.

---

## HPC Optimizations (CoreNEURON)

The NEURON engine is purpose-built for HPC environments with the following design decisions:

**Continuous Usage-Dependent Synaptic Depression**
Synaptic depression variables (`AMPAD`, `GABAAD`, `NMDAD`) are modeled as continuous ODEs solved via `cnexp`, driven by presynaptic voltage. This eliminates the numerical instabilities and queue-sorting bottlenecks caused by discrete `NET_RECEIVE` jumps in vectorized CoreNEURON builds.

**Cell-Relative Stable RNG**
Random number generation uses the `Random123` generator keyed to stable `(Cell_GID, Synapse_ID)` tuples. Results are **mathematically deterministic and bitwise identical** provided the same number of MPI ranks (`nhost`) is used each time.

**SIMD & SoA Memory Layout**
Compiled with Structure of Arrays (SoA) layout and SIMD math optimizations (`-DNRN_ENABLE_MATH_OPT=ON`) to maximize CPU cache line efficiency. In practice, this yields nearly twice the throughput of a standard NEURON build.

**Clean Terminal Output**
`sys.stdout` silencing wrappers suppress C++ build banners while displaying a real-time simulation progress bar.

---

## Key Observations

A few non-obvious behaviors encountered during development that may save time for future contributors:

**CVode is incompatible with this model**
This model drives synapses via both continuous voltage transfer *and* discrete spike events simultaneously — an uncommon pattern in NEURON. When `ParallelContext` is active, CVode fails because it cannot synchronize variable-timestep voltage transfers across MPI ranks. CVode was found to work with either continuous *or* discrete transfer alone, but not both combined. Use the fixed `dt` integrator.

**`--record-voltages` and `--coreneuron` are mutually exclusive**
CoreNEURON is designed for memory efficiency and only transfers spike data back to NEURON. Voltage recording requires the standard NEURON backend. Recompile with `nrnivmodl mod-files` (no `-coreneuron`) and run without the `--coreneuron` flag when voltage traces are needed.

**MPI scaling and optimal core allocation**
NEURON with MPI provides superior single-simulation performance. Speedup scales roughly linearly with core count until voltage/spike transfer overhead across ranks dominates — typically around 8 cores for this model size. For parameter sweeps, running each simulation on 4–8 cores and distributing jobs across the remaining cores via GNU Parallel is more efficient than saturating all cores on a single run. For example, on a 48-core machine: 11 simultaneous simulations × 4 cores each, with 4 cores reserved for OS buffer activity. Brian2 supports parallelization only via OpenMP, and scaling is poor with increasing thread count for individual simulations. For Brian2 parameter sweeps, omit OpenMP and use GNU Parallel to distribute single-threaded jobs across all physical cores.

> GPU parallelization has not been tested for this model, but both NEURON and Brian2 support it, and it is left for experimentation.

**Higher IN cell spike rate in Brian2**
With identical configuration, random seeds, and the `euler` solver, IN cells in Brian2 exhibit approximately twice the spike rate seen in NEURON. No numerical discrepancy source was identified, and the difference is likely attributable to architectural differences between the two simulators. If a closer match is required, apply an empirical scaling constant to the conductance of the AMPA synapses from `PYso→IN`.

---

## Key Issues in the Original DynaSim Code

Most published results appear correct; however, several errors in the original code and minor discrepancies with the paper were identified during replication:

**Critical bug in `netcon_nearest_neighbors`**
The original DynaSim function uses an `and` condition to compare pre- and postsynaptic population sizes in the initial `if` block. This produces an all-to-all connectivity matrix whenever the two populations differ in size. Changing the condition to `or` restores the nearest-neighbor algorithm described in the paper and is required to replicate the reported results.

**Fixed axial conductance between pyramidal compartments**
The original DynaSim code uses a fixed soma–dendrite axial conductance of 1.75 µS for all pyramidal cells. With this value (and the corrected connectivity), the propofol condition produces only an I-state with no C-state. The source paper by Benita et al. (2012), from which the cortical model was adapted, specifies that axial conductance varies normally across cells (σ = ±0.1 µS). Implementing this cell-to-cell variability restores the expected alternating C-state/I-state dynamics. The compartmental coupling equations in the supplementary material of Soplata et al. are also slightly incorrect.

**Discrepancy in Table 1 Poisson conductance**
The paper states that Poisson input conductance covaries with the AMPA conductances (`PYso→PYdr` and `TC→PYdr`) as propofol dose increases. Implementing this logic produced results inconsistent with the paper's hypothesis (C-state percentage increased with dose instead of decreasing). Table 2, by contrast, uses a fixed Poisson conductance of 0.004 mS/cm². Applying a fixed value of 0.005 mS/cm² (the original propofol-state value) to both tables recovers the expected behavior: C-state proportion decreases monotonically with increasing dose. This suggests the AMPA conductance change alone is sufficient to model propofol dose escalation. Both versions of the data are included for reference: `tables_propofol_BAD` uses the method described in the paper; `tables_propofol` uses the corrected method.

---

## Credits

- Original paper: [Rapid thalamocortical network switching mediated by cortical synchronization underlies propofol-induced EEG signatures: a biophysical model](https://doi.org/10.1152/jn.00068.2022)
- Original DynaSim code: [soplata-2023-thalcort-code](https://github.com/asoplata/soplata-2023-thalcort-code)
- The model was translated from DynaSim to NEURON and Brian2 with new analysis scripts. Code refactoring and documentation were assisted by **Gemini Pro** along with some parts of this README ;)
- Complete model with data: [thalamocortical-propofol-soplata23-drive](https://drive.google.com/file/d/1dj6pCTzdF7-JF-nQTx-jgqLy_QSZEy24/view?usp=sharing)

---

## Author

**Virendra Patel**
Master of Science by Research (MSR), Electrical Engineering
Indian Institute of Technology (IIT) Delhi
*Focus: Computational Neuroscience, Biophysical Thalamocortical Modeling, Signal Processing, Machine Learning*

Feel free to open an issue or reach out directly if you encounter problems during installation or with the model itself.

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
