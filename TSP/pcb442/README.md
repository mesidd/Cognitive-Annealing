## TRAVELING SALESMAN PROBLEM: PCB442

### KNOWN OPTIMUM: 50778
### BEST RESULT: 54634.76
### BEST TIME: 1.35 SECONDS
### NEAR-OPTIMALITY: 92.40% (7.60% gap)

<img width="1578" height="1566" alt="pcb442" src="https://github.com/user-attachments/assets/25febbb7-1c47-474a-9a59-915eb90c4ef9" />

### Terminal Output:

<img width="1178" height="342" alt="terminal_pcb" src="https://github.com/user-attachments/assets/7122d794-2965-4d76-9ac4-e68c5c374b23" />

### Key Files:
- main.py: This module contains the 316-line solver logic for the TSP application.
- pcb442.tsp: The standard TSPLIB input file used for testing.

### QUICK START (Testing the Solver)
To replicate the results, perform the following steps:
```bash
# 1. Clone the repository

git clone https://github.com/mesidd/Cognitive-Annealing.git
cd Cognitive-Annealing/TSP/pcb442

# 2. Install dependencies (Requires Python 3.x and matplotlib)

pip install matplotlib

# 3. Run the solver (Results will vary slightly due to stochastic process)

python main.py
```

### NOTE: OUR SOLVER IS A STOCHASTIC HEURISTIC. FOR AN OPTIMUM RESULT RUN IT A FEW TIMES.
