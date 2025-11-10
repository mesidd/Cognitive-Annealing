# Boolean Satisfiability (SAT)

This is **not** a standalone SAT solver. Instead, it is a specific application of the **Cognitive Annealing (CA)** framework, demonstrating how the general-purpose annealing engine can be “plugged into” the SAT problem domain.

The **core CA engine** remains unchanged from the one used to solve **TSP**, **Protein Folding**, and **Generative Chemistry**, as described in the paper.

---

## Benchmark Results

The solver replicates the results from the paper, showcasing its ability to:

- Solve satisfiable problems correctly  
- Find solutions in flawed *unsatisfiable* benchmarks  
- Converge on the true minimum for *contradictory* problems  

| Test Case | Problem (DIMACS CNF) | Final Cost (Unsatisfied Clauses) | Time (s) |
|------------|----------------------|----------------------------------|-----------|
| 1. Satisfiable | uf20-01.cnf | 0 | 0.05 s |
| 2. Flawed UNSAT | uuf50-01.cnf | 0 | 0.05 s |
| 3. True UNSAT | Built-in Example (`[[1], [-1]]`) | 1 | 0.11 s |

> *Note: Times are approximate and based on runs on standard commodity hardware.*

---

## 🧩 Framework Implementation

Following the **Cognitive Annealing** architecture, this solver is built from two domain-specific components that are “plugged into” the main, problem-agnostic engine:

### State & Mutator ("Creator")
- **State:** A dictionary of boolean (`True`/`False`) assignments for all variables.  
- **Mutator:** Randomly flips the value of a single variable to generate a neighboring state.

### Objective Function ("Metric")
- **Metric:** Counts the total number of unsatisfied clauses for the current state.  
- **Goal:** Minimize this cost to **0**.

### Terminal Output:

<img width="2266" height="1344" alt="Screenshot 2025-11-10 at 1 24 24 PM" src="https://github.com/user-attachments/assets/c0a7e3cb-ce9d-4440-84fa-018836dd1387" />

---
### QUICK START (Testing the Solver)
To replicate the results, perform the following steps:
```bash
# 1. Clone the repository

git clone [https://github.com/mesidd/Cognitive-Annealing.git](https://github.com/mesidd/Cognitive-Annealing.git)
cd Cognitive-Annealing/TSP/SAT

# 2. Run the solver (Results will vary slightly due to stochastic process)

python main.py

```

### NOTE: Our solver is a stochastic heuristic. For an optimum result, run it a few times.
