## TRAVELING SALESMAN PROBLEM: RL1889

### KNOWN OPTIMUM: 316536
### BEST RESULT: 353128.77
### BEST TIME: 12.88 SECONDS
### NEAR-OPTIMALITY: 88.44% (11.56% gap)

<img width="1568" height="1536" alt="rl1889" src="https://github.com/user-attachments/assets/45802292-0883-443c-bc95-6039bad2b58c" />

### Terminal Output:

<img width="2290" height="648" alt="terminal_rl" src="https://github.com/user-attachments/assets/4570287c-4517-4d0c-87d4-6fed174fcaef" />

### Key Files:
- main.py: This module contains the 317-line solver logic for the TSP application.
- pcb442.tsp: The standard TSPLIB input file used for testing.

### QUICK START (Testing the Solver)
To replicate the results, perform the following steps:
```bash
# 1. Clone the repository

git clone https://github.com/mesidd/Cognitive-Annealing.git
cd Cognitive-Annealing/TSP/rl1889

# 2. Install dependencies (Requires Python 3.x and matplotlib)

pip install matplotlib

# 3. Run the solver (Results will vary slightly due to stochastic process)

python main.py
```

### NOTE: OUR SOLVER IS A STOCHASTIC HEURISTIC. FOR AN OPTIMUM RESULT RUN IT A FEW TIMES.
