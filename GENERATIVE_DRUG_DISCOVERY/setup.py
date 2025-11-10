import os
import subprocess
import requests

def download_pdb(pdb_id, output_path):
    """Downloads a PDB file from the RCSB database."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_path, 'w') as f:
            f.write(response.text)
        print(f"Successfully downloaded '{pdb_id}.pdb'")
        return True
    else:
        print(f"Error: Failed to download PDB ID {pdb_id}")
        return False

def prepare_protein(pdb_file, output_pdbqt_file):
    """
    Prepares the protein target for docking using Open Babel.
    This step removes water, adds hydrogens, and converts to PDBQT format.
    NOTE: This requires Open Babel to be installed.
    """
    print(f"Preparing protein '{pdb_file}' for docking...")
    # Command to run Open Babel
    # -i pdb: input format is PDB
    # -o pdbqt: output format is PDBQT
    # -d: delete water molecules
    # -p: add hydrogens assuming neutral pH (7.4)
    command = [
        'obabel', pdb_file, 
        '-O', output_pdbqt_file, 
        '-d', 
        '-p', '7.4'
    ]
    
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        print(f"   ...Successfully created '{output_pdbqt_file}'")
    except FileNotFoundError:
        print("\n--- ERROR ---")
        print("Open Babel (`obabel`) command not found.")
        print("Please install it. For macOS with Homebrew: 'brew install open-babel'")
        print("For Linux: 'sudo apt-get install openbabel'")
        exit()
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running Open Babel: {e.stderr}")
        exit()


if __name__ == "__main__":
    print("--- Generative Drug Design: Phase 1 Setup ---")
    
    # The target protein for our drug discovery project
    TARGET_PDB_ID = '6LU7'
    TARGET_PDB_FILE = f"{TARGET_PDB_ID}.pdb"
    TARGET_PDBQT_FILE = "target_protein.pdbqt"

    # Download the PDB file if it doesn't exist
    if not os.path.exists(TARGET_PDB_FILE):
        download_pdb(TARGET_PDB_ID, TARGET_PDB_FILE)
    else:
        print(f"'{TARGET_PDB_FILE}' already exists.")

    # This requires Open Babel, a command-line tool for cheminformatics
    if os.path.exists(TARGET_PDB_FILE):
        prepare_protein(TARGET_PDB_FILE, TARGET_PDBQT_FILE)

        print("\n--- Setup Complete ---")
        print("You now have the prepared protein target 'target_protein.pdbqt'.")
        print("The next step will be to build the Creative Engine to invent molecules.")
