# Cognitive Annealing: A General Heuristic Framework

This repository contains the official source code for the paper:  
**_Cognitive Annealing: A General Heuristic Framework for Diverse NP-Hard Optimization._**

This work presents a lightweight, general-purpose heuristic framework based on **Simulated Annealing**.  
It is capable of solving a diverse set of **NP-hard optimization problems** with minimal modification.

---

## Overview

We demonstrate its power by applying the **same core engine** to four distinct domains:

- **Logistics:** Traveling Salesman Problem (TSP)  
- **Core Computation:** Boolean Satisfiability (SAT)  
- **Biophysics:** *Ab initio* Protein Folding  
- **Cheminformatics:** Generative Drug Discovery  

📄 **[Read the Full Paper on arXiv (To Be Updated)]**

---

## Core Concept: The Modular Framework

The framework's power comes from **decoupling the solver from the problem**.  
It consists of three modular components:

### 1. The Annealing Core (The Engine)
- Manages the temperature, cooling schedule, and acceptance probability.
- Universal and problem-agnostic.
- Implements the **Simulated Annealing** logic that drives all applications.

### 2. The State & Mutator (The "Creator")
- Defines what a “solution” looks like and how to create a “neighbor” solution.
- Domain-specific: Each problem provides its own version.
- Enables the solver to explore new configurations while maintaining feasibility.

### 3. The Objective Function (The "Metric")
- Scores the “quality” or “energy” of a given state.
- Determines whether a new state is better (lower energy or higher fitness).

> **To solve a new problem**, you only need to define:  
> - A new **Mutator** (how to generate a new candidate solution).  
> - A new **Metric** (how to evaluate solution quality).

---

## Case Studies & Solvers

This repository includes code for **four domains** demonstrating the flexibility of the Cognitive Annealing engine.

---

### 1️⃣ Traveling Salesman Problem (TSP)
**Description:** Solves the classic routing problem.

- **State:** An ordered list of cities.  
- **Mutator:** A 2-opt swap (reversing a sub-sequence of the tour).  
- **Metric:** Total Euclidean distance of the route.  
- **Benchmark:** `berlin52` (from TSPLIB).

---

### 2️⃣ Boolean Satisfiability (SAT)
**Description:** Solves logical constraint problems.

- **State:** A dictionary of True/False assignments for all variables.  
- **Mutator:** Flips the value of one random variable.  
- **Metric:** Number of unsatisfied clauses.  
- **Benchmark:** `uf20-01.cnf` (from SATLIB/DIMACS).

---

### 3️⃣ Ab Initio Protein Folding
**Description:** Finds the lowest-energy (most stable) 3D structure of a protein from an unfolded state.

- **State:** 3D coordinates of all atoms (managed by OpenMM).  
- **Mutator:** A short, high-temperature molecular dynamics “jiggle” (via OpenMM).  
- **Metric:** Potential energy (kJ/mol) calculated using the AMBER force field.  
- **Benchmark:** Chignolin (`PDB: 1L2Y`).

### 4️⃣ Generative Chemistry
**Description:** Discovers novel, valid molecules that are structurally similar to a target drug.

- **State:** A SMILES string (text representation of a molecule).  
- **Mutator:** Valid chemical transformations (e.g., adding rings, extending chains) using RDKit.  
- **Metric:** Tanimoto Similarity to a target molecule.  
- **Benchmark:** Target drug – *Celecoxib.*

---

## How to Apply the Framework to Your Problem

To solve a **new NP-hard problem**, define just two functions:

### 1. The State & Mutator
- **State:** What a “solution” looks like (e.g., a list of cities).  
- **Mutator:** How to make a small, random change (e.g., swap two cities).  

### 2. The Objective Function
- A metric that scores solution quality (e.g., total route distance).  
- The solver minimizes (or maximizes) this metric through iterative updates.

Once defined, plug these into the **Cognitive Annealing engine**, and the solver will automatically optimize your problem.

---

## Installation:

Note: Each case study may require additional dependencies or environment setup.
Please refer to the respective folder’s documentation for detailed installation instructions.

## 🧾 Citing This Work

If you use this framework or code in your research, please cite our paper:

```bibtex
@misc{sharma2025cognitive,
  title        = {Cognitive Annealing: A General Heuristic Framework for Diverse NP-Hard Optimization},
  author       = {Siddhartha Sharma},
  year         = {2025},
  eprint       = {YOUR_ARXIV_ID_HERE},  % e.g., 2511.12345
  archivePrefix= {arXiv},
  primaryClass = {cs.AI}
}

## 📄 License

This project is licensed under the **MIT License**.  
See the [LICENSE](./LICENSE) file for full details.
