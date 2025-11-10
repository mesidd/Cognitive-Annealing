import numpy as np
import time
import math
import random
import os
import shutil
from openmm.app import *
from openmm import *
import openmm.unit as unit
from Bio.PDB import PDBParser, Superimposer
from pdbfixer import PDBFixer

# --- Definitive Monte Carlo Folding Engine ---

class FoldingSolver:
    """A self-contained class for Monte Carlo protein folding using an MD-based mover."""

    def __init__(self, pdb_filename):
        try:
            self.pdb = PDBFile(pdb_filename)
            self.topology = self.pdb.topology
            self.forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
            self.system = self.forcefield.createSystem(self.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
            self.integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
            self.simulation = Simulation(self.topology, self.system, self.integrator)
            self.simulation.context.setPositions(self.pdb.positions)
            self.is_initialized = True
            print("1. Initializing the physics engine (OpenMM)... Success.")
        except Exception as e:
            print(f"   Error initializing solver: {e}")
            self.is_initialized = False

    def get_energy(self):
        state = self.simulation.context.getState(getEnergy=True)
        return state.getPotentialEnergy()

    def _mutate(self, temperature, num_steps):
        original_temp = self.integrator.getTemperature()
        self.integrator.setTemperature(temperature * unit.kelvin)
        self.simulation.step(num_steps)
        self.integrator.setTemperature(original_temp)

    def solve(self, params):
        if not self.is_initialized: return None

        best_pdb_file = "best_structure.pdb"
        snapshot_dir = "folding_snapshots"
        if params['snapshot_interval'] > 0:
            if os.path.exists(snapshot_dir):
                shutil.rmtree(snapshot_dir)
            os.makedirs(snapshot_dir)

        current_energy_obj = self.get_energy()
        current_energy = current_energy_obj.value_in_unit(unit.kilojoules_per_mole)
        best_energy = current_energy
        
        temperature_mc = params['initial_temp_mc']
        
        print(f"\n2. Starting Annealing. Initial Energy: {current_energy:.2f} kJ/mol")

        for i in range(params['iterations']):
            original_coords = self.simulation.context.getState(getPositions=True).getPositions()
            original_energy = current_energy
            
            self._mutate(temperature=params['move_temp'], num_steps=params['move_steps'])
            new_energy_obj = self.get_energy()
            new_energy = new_energy_obj.value_in_unit(unit.kilojoules_per_mole)
            
            k_B = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
            accept = False
            if new_energy < current_energy:
                accept = True
            elif temperature_mc > 0 and not (np.isnan(new_energy) or np.isinf(new_energy)):
                try:
                    prob = math.exp((current_energy_obj - new_energy_obj) / (k_B * temperature_mc * unit.kelvin))
                    if random.random() < prob:
                        accept = True
                except OverflowError: pass
            
            if accept:
                current_energy = new_energy
            else:
                self.simulation.context.setPositions(original_coords)
                current_energy = original_energy
            
            if current_energy < best_energy:
                best_energy = current_energy
                best_coords = self.simulation.context.getState(getPositions=True).getPositions()
                with open(best_pdb_file, 'w') as f:
                    PDBFile.writeFile(self.topology, best_coords, f)
            
            # It saves the CURRENT state, showing the entire trajectory.
            if params['snapshot_interval'] > 0 and (i+1) % params['snapshot_interval'] == 0:
                snapshot_file = os.path.join(snapshot_dir, f"snapshot_{i+1:06d}.pdb")
                current_snapshot_coords = self.simulation.context.getState(getPositions=True).getPositions()
                with open(snapshot_file, 'w') as f:
                    PDBFile.writeFile(self.topology, current_snapshot_coords, f)

            temperature_mc *= params['cooling_rate']

            if (i + 1) % 2500 == 0:
                print(f"Iter {(i+1):6d} | MC Temp: {temperature_mc:8.2f}K | Current E: {current_energy:12.2f} | Best E: {best_energy:12.2f}")
        
        print(f"\nAnnealing Complete. Best energy found: {best_energy:.2f} kJ/mol")
        print(f"Best structure saved to '{best_pdb_file}'")
        return best_pdb_file

def verify_final_result(ground_truth_raw_pdb, final_folded_pdb):
    print("\n--- Phase 5: Final Verification ---")
    clean_ref_pdb = "clean_ground_truth.pdb"
    print(f"1. Creating a clean reference from '{ground_truth_raw_pdb}'...")
    try:
        fixer = PDBFixer(filename=ground_truth_raw_pdb); fixer.findMissingResidues(); fixer.findNonstandardResidues(); fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(True); fixer.findMissingAtoms(); fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)
        with open(clean_ref_pdb, 'w') as f: PDBFile.writeFile(fixer.topology, fixer.positions, f)
        print(f"   ...Saved clean reference to '{clean_ref_pdb}'")
    except Exception as e: print(f"   Could not create clean reference PDB: {e}"); return
    print("\n2. Calculating final energies...")
    gt_solver = FoldingSolver(clean_ref_pdb)
    if not gt_solver.is_initialized: return
    ground_truth_energy = gt_solver.get_energy().value_in_unit(unit.kilojoules_per_mole)
    folded_solver = FoldingSolver(final_folded_pdb)
    if not folded_solver.is_initialized: return
    final_energy = folded_solver.get_energy().value_in_unit(unit.kilojoules_per_mole)
    print(f"   Ground Truth Energy ('{clean_ref_pdb}'): {ground_truth_energy:.2f} kJ/mol")
    print(f"   Our Best Folded Energy ('{final_folded_pdb}'): {final_energy:.2f} kJ/mol")
    print("\n3. Calculating Aligned C-alpha RMSD...")
    try:
        parser = PDBParser(QUIET=True); struct_ref = parser.get_structure('ref', clean_ref_pdb); struct_mobile = parser.get_structure('mobile', final_folded_pdb)
        ref_atoms = [atom for atom in struct_ref.get_atoms() if atom.name == 'CA']; mobile_atoms = [atom for atom in struct_mobile.get_atoms() if atom.name == 'CA']
        super_imposer = Superimposer(); super_imposer.set_atoms(ref_atoms, mobile_atoms); rmsd = super_imposer.rms
        print(f"\n   Final Aligned C-alpha RMSD: {rmsd:.2f} Angstroms")
        if rmsd < 4.0:
            print("\n   *******************************************************"); print("   *** PHENOMENAL SUCCESS! The fold is correct.      ***"); print("   *******************************************************")
    except Exception as e: print(f"   An error occurred during RMSD calculation: {e}")

if __name__ == "__main__":
    PDB_FILE_PATH = '1L2Y.pdb'
    INITIAL_STRUCTURE_PDB = 'initial_unfolded.pdb'
    
    if not os.path.exists('best_structure.pdb'):
        if not os.path.exists(INITIAL_STRUCTURE_PDB):
            print(f"Error: Run 'setup_initial_structure.py' first.")
        else:
            solver_params = {
                'initial_temp_mc': 350, 'cooling_rate': 0.9999, 'iterations': 10000,
                'move_temp': 600, 'move_steps': 25,
                'snapshot_interval': 2000
            }
            solver = FoldingSolver(INITIAL_STRUCTURE_PDB)
            solver.solve(solver_params)

    if os.path.exists('best_structure.pdb'):
        verify_final_result(PDB_FILE_PATH, 'best_structure.pdb')
    else:
        print("\n'best_structure.pdb' not found. Run the simulation to generate it.")
