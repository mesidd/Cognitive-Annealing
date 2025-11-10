# Protein Folding - Trp-cage mini-protein (PDB ID: 1L2Y)

## C-alpha RMSD: 1.96 Angstroms

This project demonstrates how the general-purpose **Cognitive Annealing (CA)** solver can be applied to complex biophysics problems.  

---

## Framework Implementation

The **Cognitive Annealing** framework decouples the **core solver** (the *Engine*) from the problem’s definition (the *State* and *Metric*). 

For this protein folding problem, the components are defined as follows:

---

### 1. The State & Mutator (The "Creator")

- **State:**  
  Represents the 3D coordinates of all atoms in the protein chain, managed by the **OpenMM** physics engine.

- **Mutator:**  
  A problem-specific *"jiggle"* function.  
  It generates a neighboring state by applying a short, **high-temperature molecular dynamics (MD)** simulation step.  
  This produces a new, **physically valid** conformation of the protein.

---

### 2. The Objective Function (The "Metric")

- **Metric:** The *cost* or *energy* of the current state.  
- Computed as the **Potential Energy** (in *kJ/mol*) of the structure, using the **AMBER** force field via OpenMM.  
- **Goal:** Minimize this energy to identify the most stable folded configuration.

---
<img width="2028" height="1338" alt="protein_image" src="https://github.com/user-attachments/assets/62b87e40-d339-407b-9b81-cc1eb3e8d7d5" />

### Terminal Output:

<img width="1208" height="742" alt="protein_terminal" src="https://github.com/user-attachments/assets/e38a6c93-401b-4d11-b33b-3c4c6e07b0ee" />

### QUICK START (Testing the Solver)
To replicate the results, perform the following steps:
```bash
# 1. Clone the repository

git clone [https://github.com/mesidd/Cognitive-Annealing.git](https://github.com/mesidd/Cognitive-Annealing.git)
cd Cognitive-Annealing/TSP/PROTEIN_FOLDING

# 2. Install dependencies (Requires Python 3.x)

# Create and activate a new environment
conda create -n proteinfold python=3.10
conda activate proteinfold

# Install key dependencies
conda install -c conda-forge openmm pdbfixer
pip install biopython numpy

# Note: Due to the complexity of molecular simulation dependencies, you may need to install additional system libraries or GPU drivers.
# It is strongly recommended to use Conda, which automatically handles most binary and compatibility issues for OpenMM and PDBFixer.

# 3. Run the solver (Results will vary slightly due to stochastic process)

python setup_initial_structure.py
python run_folding.py

```
### NOTE: Our solver is a stochastic heuristic. For an optimum result, run it a few times.

### 
