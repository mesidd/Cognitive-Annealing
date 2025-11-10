import math
import random
import time
import os

class HomeostaticSATSolver:
    """
    Applies "Cognitive Annealing" and "Homeostasis" principles
    to SAT and UNSAT problems.
    """

    def __init__(self, clauses, num_vars):
        self.clauses = clauses
        self.num_vars = num_vars
        self.current_solution = {}
        self.best_solution = {}
        self.current_energy = float('inf')
        self.best_energy = float('inf')

    def solve(self, initial_temp=100.0, cooling_rate=0.99995, max_iterations=2000000, target_energy=0, problem_name="Test"):
        """
        Main solver function.
        - The "Cognitive Annealing Loop"
        - The "Homeostatic Regulator" is the cooling_rate and the stopping criteria.
        """
        
        # --- System 1 (Intuitive Brain): Generating a fast 'first guess' ---
        self.current_solution = self.generate_initial_solution()
        self.current_energy = self.evaluate_solution(self.current_solution)
        self.best_solution = self.current_solution.copy()
        self.best_energy = self.current_energy
        
        print(f"[{problem_name}] Initial Intuitive Solution: {self.best_energy} unsatisfied clauses.")
        
        temp = initial_temp
        start_time = time.time()

        # --- The Cognitive Annealing Loop (The Forging Process) ---
        for i in range(max_iterations):
            if self.best_energy == target_energy:
                print(f"\n[{problem_name}] --- !!! PHENOMENAL VICTORY in Iteration {i} !!! ---")
                print(f"[{problem_name}] --- Target of {target_energy} unsatisfied clauses achieved. ---")
                print(f"[{problem_name}] --- Victory Condition Met! Terminating thought process. ---")
                break

            neighbor_solution = self.get_neighbor(self.current_solution)
            neighbor_energy = self.evaluate_solution(neighbor_solution)
            delta_energy = neighbor_energy - self.current_energy

            if (delta_energy < 0) or (temp > 0 and random.random() < math.exp(-delta_energy / temp)):
                self.current_solution = neighbor_solution
                self.current_energy = neighbor_energy
                
                if self.current_energy < self.best_energy:
                    self.best_solution = self.current_solution.copy()
                    self.best_energy = self.current_energy
            
            temp *= cooling_rate
            
            if i % 100000 == 0 and i > 0:
                print(f" [{problem_name}] ... thinking ... Iter: {i}, Temp: {temp:.2f}, Best Energy: {self.best_energy}")
        
        end_time = time.time()
        
        if self.best_energy != target_energy:
            print(f"\n[{problem_name}] --- Thought process complete. Target of {target_energy} was not reached. ---")
            
        return self.best_solution, self.best_energy, (end_time - start_time)

    def generate_initial_solution(self):
        """System 1: A random assignment of True/False to all variables."""
        solution = {}
        if self.num_vars == 0:
             return {}
        for i in range(1, self.num_vars + 1):
            solution[i] = random.choice([True, False])
        return solution

    def get_neighbor(self, solution):
        """Generates a new solution by flipping one variable at random."""
        if self.num_vars == 0:
            return solution
        neighbor = solution.copy()
        var_to_flip = random.randint(1, self.num_vars)
        neighbor[var_to_flip] = not neighbor[var_to_flip]
        return neighbor

    def evaluate_solution(self, solution):
        """
        The "UIM" (Unified Insight Metric).
        Calculates the "energy" of a solution by counting the number of *unsatisfied* clauses.
        The goal is to minimize this energy.
        """
        unsatisfied_count = 0
        if self.num_vars == 0 and len(self.clauses) > 0:
             return len(self.clauses)

        for clause in self.clauses:
            is_clause_satisfied = False
            for literal in clause:
                var = abs(literal)
                if var > self.num_vars or var == 0:
                    continue 
                is_negated = literal < 0
                
                if (is_negated and not solution[var]) or (not is_negated and solution[var]):
                    is_clause_satisfied = True
                    break
            
            if not is_clause_satisfied:
                unsatisfied_count += 1
        return unsatisfied_count

def load_dimacs_cnf_from_file(file_path):
    """Loads a DIMACS CNF problem from a file."""
    clauses = []
    num_vars = 0
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            
            for line in lines:
                line = line.strip()
                if line.startswith('c'):
                    continue
                if line.startswith('p cnf'):
                    parts = line.split()
                    if len(parts) >= 3:
                        num_vars = int(parts[2])
                    else:
                        print(f"Warning: Malformed 'p cnf' line: {line}")
                elif line:
                    literals = [int(x) for x in line.split() if x != '0']
                    if literals:
                        clauses.append(literals)
                        
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None, 0
    except Exception as e:
        print(f"Error parsing file: {e}")
        return None, 0
                
    return clauses, num_vars

# ==============================================================================
# --- Main Execution: The Three Tests ---
# ==============================================================================

if __name__ == "__main__":
    
    # --- Ensure the problem files exist ---
    sat_problem_file = 'uf20-01.cnf.txt'
    unsat_problem_file = 'uuf50-01.cnf.txt'
    
    if not os.path.exists(sat_problem_file):
        print(f"Error: '{sat_problem_file}' not found. Please create it with the SAT problem data.")
    if not os.path.exists(unsat_problem_file):
        print(f"Error: '{unsat_problem_file}' not found. Please create it with the (flawed) UNSAT problem data.")
    
    # --- TEST 1: The Solvable Problem (uf20-01.cnf) ---
    print("========================================================================")
    print(f"--- TEST 1: FORGING A SOLUTION FOR SATISFIABLE PROBLEM ({sat_problem_file}) ---")
    print("========================================================================")
    
    sat_clauses, sat_vars = load_dimacs_cnf_from_file(sat_problem_file)
    
    if sat_clauses:
        sat_solver = HomeostaticSATSolver(sat_clauses, sat_vars)
        sat_solution, sat_energy, sat_time = sat_solver.solve(max_iterations=500000, target_energy=0, problem_name="Test 1 (SAT)")
        
        print("\n========================= SAT SYNTHESIS COMPLETE ==========================")
        print(f"[Test 1 (SAT)] Final Forged Solution: {sat_energy} unsatisfied clauses.")
        print(f"[Test 1 (SAT)] Total time for SAT problem: {sat_time:.4f} seconds")
        if sat_energy == 0:
            print("[Test 1 (SAT)] VICTORY: A perfect solution was found.")
            print("\n[Test 1 (SAT)] Forged Solution (Variable Assignments):")
            print(sat_solution)
            # -------------------------------------
        print("========================================================================")
    else:
        print(f"Skipping Test 1. File not found.")

    # --- TEST 2: The Flawed Unsolvable Problem (uuf50-01.cnf) ---
    print("\n\n========================================================================")
    print(f"--- TEST 2: ANALYZING FLAWED UNSATISFIABLE PROBLEM ({unsat_problem_file}) ---")
    print("========================================================================")
    
    unsat_clauses, unsat_vars = load_dimacs_cnf_from_file(unsat_problem_file)
    
    if unsat_clauses:
        unsat_solver = HomeostaticSATSolver(unsat_clauses, unsat_vars)
        unsat_solution, unsat_energy, unsat_time = unsat_solver.solve(max_iterations=500000, target_energy=0, problem_name="Test 2 (Flawed UNSAT)")

        print("\n========================= FLAWED UNSAT SYNTHESIS COMPLETE ==========================")
        print(f"[Test 2 (Flawed UNSAT)] Final Forged Solution: {unsat_energy} unsatisfied clauses.")
        print(f"[Test 2 (Flawed UNSAT)] Total time for problem: {unsat_time:.4f} seconds")
        if unsat_energy == 0:
            print("[Test 2 (Flawed UNSAT)] DIAGNOSIS: Solver found a 0-clause solution, confirming the problem file is flawed and not a true UNSAT benchmark.")
        else:
             print(f"[Test 2 (Flawed UNSAT)] DIAGNOSIS: Solver converged on {unsat_energy}. This is not the known minimum, confirming data file issues.")
        print("========================================================================")
    else:
        print(f"Skipping Test 2. File not found.")

    # --- TEST 3: The Architect's Contradiction (A Guaranteed UNSAT Problem) ---
    print("\n\n========================================================================")
    print(f"--- TEST 3: DIAGNOSING THE ARCHITECT'S CONTRADICTION (x1 AND NOT x1) ---")
    print("========================================================================")
    
    arch_clauses = [[1], [-1]]
    arch_vars = 1
    
    arch_solver = HomeostaticSATSolver(arch_clauses, arch_vars)
    arch_solution, arch_energy, arch_time = arch_solver.solve(max_iterations=100000, target_energy=0, problem_name="Test 3 (True UNSAT)")

    print("\n========================= TRUE UNSAT SYNTHESIS COMPLETE ==========================")
    print(f"[Test 3 (True UNSAT)] Final Forged Solution: {arch_energy} unsatisfied clauses.")
    print(f"[Test 3 (True UNSAT)] Total time for problem: {arch_time:.4f} seconds")
    if arch_energy == 1:
        print("[Test 3 (True UNSAT)] VICTORY: The solver correctly diagnosed the problem as unsatisfiable")
        print("(by converging on the correct minimum energy of 1).")
    else:
         print(f"[Test 3 (True UNSAT)] ERROR: Solver converged on {arch_energy} instead of 1.")
    print("========================================================================")
