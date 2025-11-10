import numpy as np
import time
import math
import random
import os
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, DataStructs
from rdkit.Chem import Recap
from rdkit import rdBase

rdBase.DisableLog('rdApp.warning')

def get_starting_fragment(target_smiles):
    """Decomposes the target molecule and returns the largest fragment as a starting point."""
    print("1. Decomposing target molecule to find a smart starting point...")
    mol = Chem.MolFromSmiles(target_smiles)
    if not mol: return 'c1ccccc1'
    decomp = Recap.RecapDecompose(mol)
    fragments = list(decomp.GetAllChildren().keys())
    if not fragments:
        print("   ...No fragments found, starting with Benzene.")
        return 'c1ccccc1'
    largest_fragment = max(fragments, key=len)
    print(f"   ...Found {len(fragments)} fragments. Selecting the largest: {largest_fragment}")
    return largest_fragment

class MoleculeMutator:
    @staticmethod
    def mutate(smiles_string):
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol: return smiles_string
        operators = [MoleculeMutator._add_ring, MoleculeMutator._add_hydroxyl, MoleculeMutator._add_methyl, MoleculeMutator._add_amine, MoleculeMutator._extend_chain]
        for _ in range(10):
            operator = random.choice(operators); new_mol = operator(Chem.Mol(mol))
            if new_mol:
                try:
                    Chem.SanitizeMol(new_mol); new_smiles = Chem.MolToSmiles(new_mol)
                    if new_smiles and new_smiles != smiles_string: return new_smiles
                except: continue
        return smiles_string
    @staticmethod
    def _add_atom_to_aromatic(mol, atom_num):
        editable_mol = Chem.RWMol(mol); atoms = [a for a in editable_mol.GetAtoms() if a.GetIsAromatic() and a.GetNumExplicitHs() > 0]
        if not atoms: return None
        atom_to_modify = random.choice(atoms); new_atom_idx = editable_mol.AddAtom(Chem.Atom(atom_num)); editable_mol.AddBond(atom_to_modify.GetIdx(), new_atom_idx, Chem.BondType.SINGLE)
        return editable_mol.GetMol()
    @staticmethod
    def _add_hydroxyl(mol): return MoleculeMutator._add_atom_to_aromatic(mol, 8)
    @staticmethod
    def _add_methyl(mol): return MoleculeMutator._add_atom_to_aromatic(mol, 6)
    @staticmethod
    def _add_amine(mol): return MoleculeMutator._add_atom_to_aromatic(mol, 7)
    @staticmethod
    def _extend_chain(mol):
        editable_mol = Chem.RWMol(mol); atoms = [a for a in editable_mol.GetAtoms() if not a.GetIsAromatic() and a.GetDegree() < 4 and a.GetAtomicNum() == 6]
        if not atoms: return None
        atom_to_modify = random.choice(atoms); new_atom_idx = editable_mol.AddAtom(Chem.Atom(6)); editable_mol.AddBond(atom_to_modify.GetIdx(), new_atom_idx, Chem.BondType.SINGLE)
        return editable_mol.GetMol()
    @staticmethod
    def _add_ring(mol):
        editable_mol = Chem.RWMol(mol); atoms = [a for a in editable_mol.GetAtoms() if a.GetDegree() > 1 and a.GetTotalNumHs() > 0]
        if len(atoms) < 2: return None
        a1, a2 = random.sample(atoms, 2)
        if editable_mol.GetBondBetweenAtoms(a1.GetIdx(), a2.GetIdx()) is not None: return None
        path = Chem.GetShortestPath(mol, a1.GetIdx(), a2.GetIdx())
        if len(path) > 0 and len(path) < 4: return None
        editable_mol.AddBond(a1.GetIdx(), a2.GetIdx(), Chem.BondType.SINGLE)
        return editable_mol.GetMol()

class SimilarityEvaluator:
    def __init__(self, target_smiles):
        self.target_mol = Chem.MolFromSmiles(target_smiles);
        if not self.target_mol: raise ValueError("Invalid target SMILES string.")
        self.fpgen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
        self.target_fp = self.fpgen.GetFingerprint(self.target_mol)
    def evaluate(self, candidate_smiles):
        candidate_mol = Chem.MolFromSmiles(candidate_smiles)
        if not candidate_mol: return 0.0
        candidate_fp = self.fpgen.GetFingerprint(candidate_mol)
        return DataStructs.TanimotoSimilarity(self.target_fp, candidate_fp)

class HomeostaticSolver:
    def solve(self, params):
        print("\n2. Initializing Homeostatic Solver...")
        evaluator = SimilarityEvaluator(target_smiles=params['target_smiles'])
        current_smiles = params['start_smiles']
        current_score = evaluator.evaluate(current_smiles)
        best_smiles, best_score = current_smiles, current_score
        temperature = params['initial_temp']; stagnation_counter = 0
        print(f"   Starting with '{current_smiles}'. Initial Score: {current_score:.4f}")
        for i in range(params['iterations']):
            new_smiles = MoleculeMutator.mutate(current_smiles)
            new_score = evaluator.evaluate(new_smiles)
            if new_score > current_score or (temperature > 0 and random.random() < math.exp((new_score - current_score) / temperature)):
                current_smiles, current_score = new_smiles, new_score
            if current_score > best_score:
                best_score, best_smiles = current_score, current_smiles
                stagnation_counter = 0
            else:
                stagnation_counter += 1
            if stagnation_counter >= params['reheat_threshold']:
                print(f"--- Stagnation detected at iter {i+1}! Reheating. ---")
                temperature = params['initial_temp'] * 0.5
                stagnation_counter = 0
            temperature *= params['cooling_rate']
            if (i + 1) % 1000 == 0:
                print(f"Iter {(i+1):6d} | Temp: {temperature:.4f} | Current Score: {current_score:.4f} | Best Score: {best_score:.4f}")
        print(f"\nAnnealing Complete. Best score found: {best_score:.4f}")
        print(f"Best molecule invented: {best_smiles}")
        return best_smiles, best_score

if __name__ == "__main__":
    print(f"--- Generative Drug Design: Upgraded Homeostatic Solver ---")
    TARGET_DRUG_SMILES = "CS(=O)(=O)C1=CC=C(C=C1)C2=C(C(=O)N2C)C3=CC=C(C=C3)OC(F)(F)F"
    params = { 'target_smiles': TARGET_DRUG_SMILES, 
      'start_smiles': get_starting_fragment(TARGET_DRUG_SMILES),      
      'initial_temp': 0.2, 
      'cooling_rate': 0.9, 
      'iterations': 4000, 
      'reheat_threshold': 1000 
    }
    solver = HomeostaticSolver()
    best_molecule_smiles, best_molecule_score = solver.solve(params)

    print("\n--- Final Analysis ---")
    print(f"Target Molecule: {params['target_smiles']}")
    print(f"AI Invented Molecule: {best_molecule_smiles}")
    print(f"Final Tanimoto Similarity: {best_molecule_score:.4f}")

    if best_molecule_score > 0.7: print("\n*** PHENOMENAL SUCCESS! The AI has invented a novel and highly similar molecule. ***")
    elif best_molecule_score > 0.5: print("\n*** Excellent result! The AI has discovered a promising molecular scaffold. ***")

    print("\nGenerating visualization...")

    target_mol = Chem.MolFromSmiles(params['target_smiles']); invented_mol = Chem.MolFromSmiles(best_molecule_smiles)
    img = Draw.MolsToGridImage([target_mol, invented_mol], legends=['Target Drug', 'AI Invented Molecule'], subImgSize=(400, 400))

    img.save('drug_discovery_result.png'); print("...Final comparison image saved to 'drug_discovery_result.png'")
