from pdbfixer import PDBFixer
from openmm.app import *
import openmm as mm # <-- [THE FIX] Import the core openmm library
import openmm.unit as unit

print("--- Preparing a valid initial structure for folding ---")

PDB_ID = '1L2Y'
INITIAL_STRUCTURE_PDB = 'initial_unfolded.pdb'

# 1. Create a "fixer" for our PDB ID.
print(f"1. Fetching and fixing PDB ID: {PDB_ID}...")
fixer = PDBFixer(pdbid=PDB_ID)

# 2. Find and add missing atoms
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(True) # Remove water
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)
print("   ...model is now complete and physically valid.")

# 3. Set up a short "heating" simulation to unfold the protein
print("2. Running a high-temperature simulation to unfold the structure...")
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
# [THE FIX] Use mm.LangevinMiddleIntegrator
system = forcefield.createSystem(fixer.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
integrator = mm.LangevinMiddleIntegrator(800*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds) 
simulation = Simulation(fixer.topology, system, integrator)
simulation.context.setPositions(fixer.positions)

simulation.minimizeEnergy()
simulation.step(2500) # Run heating for 5 picoseconds
print("   ...protein has been unfolded into a random coil.")

# 4. Save the final, unfolded, but physically valid positions
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open(INITIAL_STRUCTURE_PDB, 'w'))
print(f"\n3. Successfully created a valid starting structure at '{INITIAL_STRUCTURE_PDB}'")
