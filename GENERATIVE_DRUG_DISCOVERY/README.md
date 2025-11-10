# Generative Drug Discovery with Cognitive Annealing

## Targeting SARS-CoV-2 Main Protease (PDB: 6LU7)

>Tanimoto Similarity: 0.628. Interpretation: The molecule shares approximately 63% similarity with the reference ligand, suggesting comparable structural features.


This code demonstrates how the **general-purpose Cognitive Annealing framework** can be applied to a **de novo drug discovery problem**.  

The goal is to navigate a vast chemical design space to *invent* novel, valid molecules that are structurally similar to a known target drug.

---

## How It Works: The Framework

This is **not a solver built from scratch for chemistry**. Instead, it is an **application of the general Cognitive Annealing engine**. We simply *plug in* two problem-specific components:

---

### 1. The State & Mutator (The “Creator”)

- **State:**  
  A potential solution is represented by a **SMILES string**, a text-based representation of a molecule.

- **Mutator:**  
  A function that uses the **RDKit** library to apply valid chemical transformations —  
  e.g. adding a ring, extending a chain, or modifying an atom —  
  to the current SMILES string, thereby creating a new, valid “neighbor” molecule.

---

### 2. The Objective Function (The “Metric”)

- **Metric:**  
  The **Tanimoto Similarity Score**.

- **Goal:**  
  The engine aims to **maximize** this score, which measures the *structural similarity*  
  between a generated molecule and a target drug.

---

<img width="1200" height="400" alt="drug_discovery_result" src="https://github.com/user-attachments/assets/0efd6715-6b09-4448-8716-5f12a3cb307c" />

### Terminal Output:

<img width="1346" height="594" alt="Drug Discovery" src="https://github.com/user-attachments/assets/96bc80c5-9154-4c1f-b1c0-d7b33ae31566" />

### QUICK START (Testing the Solver)
To replicate the results, perform the following steps:
```bash
# 1. Clone the repository

git clone https://github.com/mesidd/Cognitive-Annealing.git
cd Cognitive-Annealing/TSP/GENERATIVE_DRUG_DISCOVERY

# 2. Install dependencies

pip install rdkit-pypi requests

# OR

conda create -n drugdiscovery python=3.10
conda activate drugdiscovery
conda install -c conda-forge rdkit
pip install requests

# Depending on your setup, additional dependencies (e.g., visualization or cheminformatics libraries) may 
# be required for full functionality.

# 3. Run the solver (Results will vary slightly due to stochastic process)

python setup.py
python generate_drug.py

```

## Summary

This experiment showcases the **problem-agnostic design** of Cognitive Annealing —  
the same optimization engine used for **TSP**, **SAT solving**, and **Protein Folding**  
can also be applied to **molecular generation** and **drug discovery**.

### NOTE: Our solver is a stochastic heuristic. For an optimum result, run it a few times.
