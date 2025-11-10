import math
import random
import time
import re
import matplotlib.pyplot as plt
from multiprocessing import Pool, Manager

# --- Core Data Structures ---

class City:
    """Represents a city with (id, x, y) coordinates."""
    def __init__(self, id, x, y):
        self.id = int(id)
        self.x = float(x)
        self.y = float(y)

    def distance_to(self, other_city):
        """Calculates Euclidean distance to another city."""
        return math.sqrt((self.x - other_city.x)**2 + (self.y - other_city.y)**2)

    def __repr__(self):
        return f"City {self.id}"

class Tour:
    """Represents a tour (a sequence of cities) and its total distance."""
    def __init__(self, cities):
        self.tour = list(cities)
        self.distance = self.calculate_distance()

    def calculate_distance(self):
        """Calculates the total distance of the tour."""
        total_distance = 0
        for i in range(len(self.tour)):
            from_city = self.tour[i]
            to_city = self.tour[(i + 1) % len(self.tour)]
            total_distance += from_city.distance_to(to_city)
        return total_distance

# --- "Cognitive Annealing" CORE ARCHITECTURE ---

class TSPCreativeBrain:
    """
    Creates novel 'proto-ideas' (neighbor solutions) 
    through chaotic, non-obvious mutations. This embodies the "Stochastic 
    Resonance" principle.
    """
    @staticmethod
    def generate_neighbor_solution(tour):
        """Creates a new tour by applying a mutation (2-opt or swap)."""
        new_tour_cities = list(tour.tour)
        mutation_type = random.random()
        
        # 70% chance for a 2-opt swap (inversion)
        if mutation_type < 0.7:
            i, j = sorted(random.sample(range(len(new_tour_cities)), 2))
            new_tour_cities[i:j] = reversed(new_tour_cities[i:j])
        # 30% chance for a single city swap
        else:
            city_to_move = new_tour_cities.pop(random.randrange(len(new_tour_cities)))
            new_tour_cities.insert(random.randrange(len(new_tour_cities)), city_to_move)
            
        return Tour(new_tour_cities)

class TSPLogicalBrain:
    """
    Applies rigorous, structured critique to a 
    solution to find verifiable improvements. This is the "Ordered Critique" 
    that complements the chaotic creation.
    """
    @staticmethod
    def apply_2_opt_refinement(tour):
        """A key function of the "Feynman Engine": finding verifiable truth via 2-opt."""
        improved = True
        tour_cities = list(tour.tour)
        while improved:
            improved = False
            num_cities = len(tour_cities)
            for i in range(num_cities - 2):
                for j in range(i + 2, num_cities):
                    p1, p2 = tour_cities[i], tour_cities[i+1]
                    p3, p4 = tour_cities[j], tour_cities[(j+1)%num_cities]
                    
                    if (p1.distance_to(p3) + p2.distance_to(p4)) < (p1.distance_to(p2) + p3.distance_to(p4)):
                        tour_cities[i+1:j+1] = reversed(tour_cities[i+1:j+1])
                        improved = True
        return Tour(tour_cities)

class HomeostaticSolver:
    """
    The core "Cognitive Annealing Loop" and "Homeostatic Regulator."
    It guides the search from a high-energy "chaotic" state (high temp)
    to a low-energy "crystallized" state (low temp).
    """
    def __init__(self, cities, params, accuracy_target, stop_event, victory_tours):
        self.cities = cities
        self.params = params
        self.accuracy_target = accuracy_target
        self.stop_event = stop_event
        self.victory_tours = victory_tours  

    def solve_trial(self, trial_num):
        """The thought process for a single parallel trial, striving for victory."""
        # Start from a random "proto-idea" (replaces System 1)
        current_tour = Tour(random.sample(self.cities, len(self.cities)))
        best_tour = current_tour
        temperature = self.params['initial_temp']
        
        stagnation_counter = 0
        iteration = 0
        progress_update_interval = 200000

        while temperature > self.params['min_temp'] and not self.stop_event.is_set():
            iteration += 1
            
            # 1. Heat (Chaotic Creation)
            new_tour = TSPCreativeBrain.generate_neighbor_solution(current_tour)
            
            # 2. Measurement (UIM Proxy)
            if self.acceptance_probability(current_tour.distance, new_tour.distance, temperature) > random.random():
                current_tour = new_tour
            
            # 3. Update Best
            if current_tour.distance < best_tour.distance:
                best_tour = Tour(current_tour.tour)
                stagnation_counter = 0 # Reset stagnation
                
                # 4. Victory Condition (Goal Met)
                if best_tour.distance <= self.accuracy_target:
                    print(f"\n--- [TRIAL {trial_num}] Target Achieved ---")
                    print(f"--- Target: <= {self.accuracy_target:.2f} | Result: {best_tour.distance:.2f} ---")
                    self.victory_tours.append(best_tour)
                    self.stop_event.set()
                    break
            else:
                stagnation_counter += 1

            # 5. Homeostatic Regulation (Reheating)
            if stagnation_counter >= self.params['max_stagnation']:
                temperature += self.params['initial_temp'] * 0.1 # Reheat
                stagnation_counter = 0

            # 6. Cooling Schedule
            temperature *= self.params['cooling_rate']
            
            if iteration % progress_update_interval == 0:
                 print(f"  ... [TRIAL {trial_num}] Processing... Iter: {iteration}, Temp: {temperature:.2f}, Best: {best_tour.distance:.2f}")

        # Final "Ordered Critique" on the best tour found
        return TSPLogicalBrain.apply_2_opt_refinement(best_tour)

    def acceptance_probability(self, current_dist, new_dist, temperature):
        """The core of the annealing decision logic."""
        if new_dist < current_dist: return 1.0
        if temperature <= 0.0: return 0.0
        return math.exp((current_dist - new_dist) / temperature)

# --- Utilities ---

def parse_tsp_file(filename):
    """
    Parses a .tsp file (TSPLIB format) and returns a list of City objects.
    Handles scientific notation for coordinates.
    """
    cities = []
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        print(f"Please make sure '{filename}' is in the same directory as the script.")
        return None

    start_index = 0
    for i, line in enumerate(lines):
        if line.strip() == "NODE_COORD_SECTION":
            start_index = i + 1
            break
    
    if start_index == 0:
        raise ValueError("Could not find NODE_COORD_SECTION in the TSP file.")

    # Regex to capture "id x y", robust to scientific notation (e.g., 2.00000e+02)
    coord_regex = re.compile(r"^\s*(\d+)\s+([\d\.e\+\-]+)\s+([\d\.e\+\-]+)\s*$")

    for line in lines[start_index:]:
        line_strip = line.strip()
        if line_strip == "EOF":
            break
        
        match = coord_regex.match(line_strip)
        if match:
            city_id, x, y = match.groups()
            cities.append(City(city_id, x, y))
            
    if not cities:
        raise ValueError("No city coordinates were parsed. Check the file format.")
        
    print(f"--- Successfully parsed {len(cities)} cities from {filename} ---")
    return cities

def plot_tour(final_tour, title):
    """Plots a single, optimized tour."""
    fig, ax = plt.subplots(figsize=(8, 8))
    fig.suptitle(title, fontsize=16)
    
    ax.set_title(f"Final Optimized Tour\nDistance: {final_tour.distance:.2f}")
    
    # Get coordinates for plotting
    x_coords = [city.x for city in final_tour.tour]
    x_coords.append(final_tour.tour[0].x) # Add start city to end
    
    y_coords = [city.y for city in final_tour.tour]
    y_coords.append(final_tour.tour[0].y) # Add start city to end

    ax.plot(x_coords, y_coords, 'o-')
    
    # Annotate city IDs
    for city in final_tour.tour:
        ax.text(city.x, city.y, str(city.id), fontsize=8, ha='right')
        
    ax.set_xlabel("X Coordinate")
    ax.set_ylabel("Y Coordinate")
    ax.grid(True)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# --- Wrapper for Parallel Execution ---
def run_solver_trial_wrapper(args):
    """Helper function to unpack arguments for the pool map."""
    solver, trial_num = args
    return solver.solve_trial(trial_num)

# --- Main Execution ---
if __name__ == "__main__":
    
    TSP_FILE = 'pcb442.tsp'
    
    print("\n" + "="*20 + f" Initializing TSP Solver for {TSP_FILE} " + "="*20)

    cities = parse_tsp_file(TSP_FILE)
    
    # Exit gracefully if file wasn't found
    if cities is None:
        exit()
    
    known_optimum = 50778
    accuracy_goal_percentage = .96 # 96% of optimal
    accuracy_target_distance = known_optimum / accuracy_goal_percentage
    
    solver_params = {
        'initial_temp': 1000,
        'cooling_rate': 0.93, 
        'min_temp': .001,
        'max_stagnation': 400
    }
    
    num_trials = 4
    print(f"--- Starting Parallel Optimization with {num_trials} trials ---")
    print(f"--- Target: Achieve a tour distance <= {accuracy_target_distance:.2f} ---")
    
    start_time = time.time()
    
    with Manager() as manager:
        stop_event = manager.Event()
        victory_tours = manager.list()
        solver = HomeostaticSolver(cities, solver_params, accuracy_target_distance, stop_event, victory_tours)
        
        pool = Pool(processes=num_trials)
        trial_args = [(solver, i + 1) for i in range(num_trials)]
        
        async_results = [pool.apply_async(run_solver_trial_wrapper, args=[arg]) for arg in trial_args]
        
        pool.close()
        
        while not stop_event.is_set():
            if all(res.ready() for res in async_results):
                break 
            time.sleep(0.1)
        
        if stop_event.is_set():
            print("--- Target Achieved! Terminating remaining trials. ---")
            pool.terminate()

        pool.join()

        all_candidate_tours = []
        for res in async_results:
            if res.ready() and res.successful():
                all_candidate_tours.append(res.get())
        
        all_candidate_tours.extend(list(victory_tours))

        if not all_candidate_tours:
             raise Exception("No trials completed successfully and no victory tour was recorded. Check parameters.")

        best_overall_tour = min(all_candidate_tours, key=lambda tour: tour.distance)
        end_time = time.time()

    print("\n" + "="*27 + " Optimization Complete " + "="*27)
    print(f"Best Tour Distance (from all trials): {best_overall_tour.distance:.2f}")
    print(f"Known Optimum for {TSP_FILE}: {known_optimum}")
    
    # Use your accuracy calculation
    accuracy_achieved = 100 * (1 - (best_overall_tour.distance - known_optimum) / known_optimum)
    print(f"This result is {accuracy_achieved:.2f}% accurate vs. the optimum.")
    
    if best_overall_tour.distance <= accuracy_target_distance:
        print(f"SUCCESS: Final result of {best_overall_tour.distance:.2f} is within the target of <= {accuracy_target_distance:.2f}.")
    else:
        print(f"FAILURE: Final result {best_overall_tour.distance:.2f} missed the target of <= {accuracy_target_distance:.2f}.")
    
    print(f"Total Computation Time: {end_time - start_time:.2f} seconds")
    print("="*72)
    
    # Plot the single best tour
    plot_tour(best_overall_tour, f"TSP Solution: {TSP_FILE}")
