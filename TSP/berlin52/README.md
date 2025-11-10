## TRAVELING SALESMAN PROBLEM: BERLIN52

### KNOWN OPTIMUM: 7542
### BEST RESULT: 7544.37
### BEST TIME: 1.88 SECONDS
### NEAR-OPTIMALITY: 99.97% (0.03% gap)


<img width="1000" height="1000" alt="berlin52" src="https://github.com/user-attachments/assets/f3204b3f-bf6f-424a-9949-4e424c99cdb7" />

### Terminal Output:

<img width="1208" height="446" alt="terminal" src="https://github.com/user-attachments/assets/a455c1f4-9060-4973-ba42-5b38b0a11add" />

### Key Files:
- main.py: This module contains the 327-line solver logic for the TSP application.
- berlin52.tsp: The standard TSPLIB input file used for testing.

### QUICK START (Testing the Solver)
To replicate the results, perform the following steps:
```bash
# 1. Clone the repository

git clone [https://github.com/mesidd/Cognitive-Annealing.git](https://github.com/mesidd/Cognitive-Annealing.git)
cd Cognitive-Annealing/TSP/berlin52

# 2. Install dependencies (Requires Python 3.x and matplotlib)

pip install matplotlib

# 3. Run the solver (Results will vary slightly due to stochastic process)

python main.py
```

### NOTE: OUR SOLVER IS A STOCHASTIC HEURISTIC. FOR AN OPTIMUM RESULT RUN IT A FEW TIMES.
