### 📌 Interactive Receptor & Grid Box Preparation Pipeline (AutoDock Vina / SWISS-MODEL)[cite: 4]

[![View Script](https://img.shields.io/badge/VIEW_SCRIPT-PREP__DOCKING.PY-blue?style=for-the-badge&logo=python)](#)
[![Status](https://img.shields.io/badge/STATUS-EXECUTABLE-brightgreen?style=for-the-badge)](#)

<details>
<summary>🔍 <b>Click to View Script Details & Usage Header</b></summary>

```python
# ==============================================================================
# SCRIPT NAME : prep_docking_advanced.py
# AUTHOR      : K. R. Janahan
# DESCRIPTION : Interactive PDB receptor preparation script for molecular docking. 
#               Extracts chains/ligands, builds AutoDock Vina grid box configs, 
#               performs loop repairs via SWISS-MODEL API, and converts to PDBQT.
# DEPENDENCIES: Biopython, PDBFixer, OpenMM, MGLTools, Requests, Python3
# USAGE       : python3 prep_docking_advanced.py protein.pdb
# ==============================================================================
```
</details>

<details>
<summary>💻 <b>Click to View / Copy Full Source Code</b></summary>

```python

import sys
import numpy as np
import time
import requests
import json
import os
import gzip
import warnings
import subprocess
from Bio.PDB import PDBParser, PDBIO, Select
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# --- JANAHAN OS // UNIVERSAL STARTUP ANIMATION BLOCK ---
try:
    import sys, os
    # 1. Automatically inject your custom tool folder into the path
    custom_path = os.path.expanduser("~/intro_script")
    if custom_path not in sys.path:
        sys.path.insert(0, custom_path)

    # 2. Dynamically detect the name of the script running
    try:
        script_name = os.path.splitext(os.path.basename(__file__))[0].replace("_", " ").upper()
    except NameError:
        script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0].replace("_", " ").upper() if sys.argv and sys.argv[0] else "AUTOMATION PIPELINE"

    # 3. Trigger animation with dynamic naming
    from bio_intro import play_intro
    play_intro(f"PROTOCOL: {script_name} // EXECUTION MATRIX")

except Exception as e:
    # DEBUG FLAG: Print the exact reason why it fails instead of silently passing!
    print(f"\n[JANAHAN OS NOTICE] Startup animation skipped: {e}\n")
# -------------------------------------------------------

# Suppress annoying Biopython renumbering warning floods
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def main():
    if len(sys.argv) != 2:
        print("Usage: python prep_docking.py <raw_protein.pdb>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure("protein", pdb_file)
    except Exception as e:
        print(f"Error loading PDB file: {e}")
        sys.exit(1)

    # 1. Identify and select Protein Chains
    all_chains = []
    for model in structure:
        for chain in model:
            if chain.id not in all_chains:
                all_chains.append(chain.id)
                
    print("\n--- FOUND PROTEIN CHAINS ---")
    print(f"Chains present in file: {', '.join(all_chains)}")
    keep_chains_input = input("Enter the chain ID(s) you want to KEEP (comma-separated, e.g., A or A,B) or press Enter to keep ALL: ").strip()
    
    if keep_chains_input:
        kept_chains_list = [c.strip().upper() for c in keep_chains_input.split(',')]
    else:
        kept_chains_list = all_chains

    print(f"-> Proceeding with chains: {', '.join(kept_chains_list)}")

    # Dictionaries to store found non-standard residues
    ligands_and_cofactors = []
    waters = []

    # 2. Scan the structure
    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = residue.id
                if res_id[0] != ' ':
                    res_name = residue.resname.strip()
                    chain_id = chain.id
                    res_num = res_id[1]
                    
                    item_info = {
                        "name": res_name,
                        "chain": chain_id,
                        "num": res_num,
                        "residue_obj": residue
                    }
                    
                    if res_name == 'HOH' or res_id[0] == 'W':
                        waters.append(item_info)
                    else:
                        ligands_and_cofactors.append(item_info)

    # 3. Display Ligands/Cofactors/Ions
    print("\n--- FOUND LIGANDS, IONS, AND COFACTORS ---")
    if not ligands_and_cofactors:
        print("None found.")
    else:
        for i, lig in enumerate(ligands_and_cofactors):
            print(f"[{i}] {lig['name']} (Chain {lig['chain']}, Residue {lig['num']})")

    # 4. Ask user for Grid Box reference ligand
    ref_ligand_idx = input("\nEnter the number [ ] of the ligand to use as the Grid Box center: ").strip()
    ref_ligand = None
    if ref_ligand_idx.isdigit() and 0 <= int(ref_ligand_idx) < len(ligands_and_cofactors):
        ref_ligand = ligands_and_cofactors[int(ref_ligand_idx)]
    else:
        print("Invalid selection. Exiting.")
        sys.exit(1)

    # 5. Ask user which items to keep
    print("\n--- PREPARING THE RECEPTOR ---")
    keep_indices = input("Enter the numbers [ ] of any OTHER ligands/ions/cofactors to KEEP in the receptor (comma-separated, or press Enter for none): ").strip()
    
    kept_items = []
    if keep_indices:
        for idx in keep_indices.split(','):
            idx = idx.strip()
            if idx.isdigit() and 0 <= int(idx) < len(ligands_and_cofactors):
                kept_items.append(ligands_and_cofactors[int(idx)]["residue_obj"])

    print(f"\nFound {len(waters)} water molecules.")
    keep_waters = input("Do you want to keep ANY important water molecules? (y/n): ").strip().lower()
    
    if keep_waters == 'y':
        print("Please enter the specific residue numbers of the waters to keep, separated by commas.")
        water_nums = input("Water numbers to keep (or press Enter to discard all): ").strip()
        if water_nums:
            allowed_water_nums = [int(n.strip()) for n in water_nums.split(',') if n.strip().isdigit()]
            for w in waters:
                if w['num'] in allowed_water_nums:
                    kept_items.append(w['residue_obj'])

    # 6. Calculate Grid Box
    print(f"\n--- CALCULATING GRID BOX FOR {ref_ligand['name']} ---")
    coords = [atom.coord for atom in ref_ligand["residue_obj"].get_atoms()]
    coords = np.array(coords)

    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)

    # FIX 1: Use the spatial center of the bounding box, NOT the mean (average).
    # This guarantees the box is perfectly perfectly centered around the extremes of the ligand.
    center = (max_coords + min_coords) / 2.0
    
    # FIX 2: Ask for padding so you can control how tight the box is. 
    # A padding of 4.0 adds 2.0 Angstroms to each side.
    try:
        padding_input = input("Enter padding for the grid box in Angstroms (Press Enter for a tight fit of 4.0): ").strip()
        padding = float(padding_input) if padding_input else 4.0
    except ValueError:
        print("Invalid input. Using default padding of 4.0")
        padding = 4.0

    dimensions = (max_coords - min_coords) + padding
    
    # Use ceil to round up to the nearest whole integer to ensure the ligand is never clipped
    size_x = int(np.ceil(dimensions[0]))
    size_y = int(np.ceil(dimensions[1]))
    size_z = int(np.ceil(dimensions[2]))

    print(f"Center (X, Y, Z): {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}")
    print(f"Size (X, Y, Z) in Angstroms: {size_x}, {size_y}, {size_z}")

    # 7. Save the Cleaned Receptor PDB 
    class ReceptorSelect(Select):
        def __init__(self, allowed_chains, allowed_items):
            self.allowed_chains = allowed_chains
            self.allowed_items = allowed_items

        def accept_residue(self, residue):
            chain_id = residue.get_parent().id
            if residue.id[0] == ' ':
                return chain_id in self.allowed_chains
            if residue in self.allowed_items:
                return True
            return False

    cleaned_pdb = f"{pdb_file.split('.')[0]}_clean.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(cleaned_pdb, ReceptorSelect(kept_chains_list, kept_items))
    print(f"\n[SUCCESS] Cleaned receptor saved to: {cleaned_pdb}")

    # =====================================================================
    # MOVED STEP 8: Generate conf.txt for AutoDock Vina (Right after Step 7)
    # =====================================================================
    repaired_pdb = f"{pdb_file.split('.')[0]}_receptor.pdb" # Pre-define for config map
    config_filename = "conf.txt"
    with open(config_filename, "w") as f:
        f.write("# Copyrights @Insilicomics Research Pvt, Ltd\n\n")
        
        f.write(f"receptor = {repaired_pdb}qt\n\n") 
        
        f.write(f"center_x = {center[0]:.3f}\n")
        f.write(f"center_y = {center[1]:.3f}\n")
        f.write(f"center_z = {center[2]:.3f}\n\n")
        
        f.write(f"size_x = {size_x}\n")
        f.write(f"size_y = {size_y}\n")
        f.write(f"size_z = {size_z}\n\n")
        
        f.write("energy_range = 4\n\n")
        
        f.write("#receptor - Specifies the receptor file in PDBQT format\n")
        f.write("#center_  - Defines the center of the search space in 3D space (in Ångströms)\n")
        f.write("#size     - Sets the size of the search space in each dimension (in Ångströms)\n")
        f.write("#spacing  - Explicitly set to 1.000 to match AutoDockTools visualizations\n")
        
    print(f"[SUCCESS] Docking configuration saved to: {config_filename}")
    print(f"\n[NOTE for AutoDockTools]: To visualize this box correctly in the GUI,")
    print(f"ensure you change the 'Spacing (angstrom)' slider from 0.375 to 1.000.")

    # NEW: Interactive checkpoint to continue or stop pipeline execution
    proceed_input = input("\nDo you want to continue with the generation of the receptor file? (y/n): ").strip().lower()
    if proceed_input != 'y':
        print("\n[INFO] Script stopped by user request. Docking configuration preserved.")
        sys.exit(0)

    # ==========================================
    # THE MAGIC FIX: Restore the SEQRES Blueprint
    # ==========================================
    seqres_lines = []
    with open(pdb_file, 'r') as raw_f:
        # ... rest of your script (SEQRES recovery and SWISS-MODEL API execution) continues normally ...
        for line in raw_f:
            if line.startswith("SEQRES"):
                # Column 12 in a PDB file holds the Chain ID
                chain_id = line[11] 
                if chain_id in kept_chains_list:
                    seqres_lines.append(line)

    if seqres_lines:
        with open(cleaned_pdb, 'r') as clean_f:
            clean_content = clean_f.read()
        with open(cleaned_pdb, 'w') as clean_f:
            clean_f.writelines(seqres_lines)
            clean_f.write(clean_content)
            
    print("[INFO] Restored SEQRES blueprint for automatic loop building.")

    # ==========================================
    # NEW STEP: Auto-Repair with SWISS-MODEL API
    # ==========================================
    repaired_pdb = f"{pdb_file.split('.')[0]}_receptor.pdb"
    print("\n--- REPAIRING MISSING LOOPS & ATOMS VIA SWISS-MODEL ---")
    
    SM_TOKEN = "56e38e7625b71da955432de8c54a1602a0b6a1c0"  
    pdb_id = os.path.basename(pdb_file)[:4].upper()
    
    try:
        if SM_TOKEN == "56e38e7625b71da955432de8c54a1602a0b6a1c0" and "YOUR" in SM_TOKEN:
            raise ValueError("Token configuration string placeholder mismatch.")

        # --- A. Fetch the FASTA sequence from RCSB ---
        print(f"Fetching FASTA sequence for {pdb_id} from RCSB PDB...")
        fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
        fasta_response = requests.get(fasta_url)
        
        if fasta_response.status_code != 200:
            raise Exception(f"Failed to fetch FASTA for {pdb_id}. HTTP {fasta_response.status_code}")
            
        fasta_lines = fasta_response.text.strip().split('\n')
        target_sequence = "".join([line for line in fasta_lines if not line.startswith(">")])
        print(f"Successfully retrieved sequence ({len(target_sequence)} amino acids).")

        # --- B. Submit the job via User Template API ---
        print("Uploading template data and submitting to SWISS-MODEL...")
        headers = {
            "Authorization": f"Token {SM_TOKEN}",
            "Content-Type": "application/json"
        }
        
        # Read the clean PDB file contents to text format as requested by the API
        with open(cleaned_pdb, "r") as f:
            template_coordinates = f.read()
            
        # Unified JSON payload schema for User-Uploaded templates
        payload = {
            "target_sequences": [target_sequence],
            "template_coordinates": template_coordinates,
            "project_title": f"{pdb_id}_Docking_Prep"
        }
        
        response = requests.post(
            "https://swissmodel.expasy.org/user_template", 
            headers=headers, 
            json=payload
        )
        
        if response.status_code not in [200, 201, 202]:
            raise Exception(f"SWISS-MODEL API submission failed with HTTP {response.status_code}: {response.text[:200]}")
            
        project_data = response.json()
        project_id = project_data.get("project_id")
        if not project_id:
            raise Exception("Did not receive a valid project ID from the server.")
            
        print(f"Job successfully submitted! Project ID: {project_id}")
        
        # --- C. Poll the server until the job is done ---
        status_url = f"https://swissmodel.expasy.org/project/{project_id}/models/summary/"
        is_finished = False
        
        print("Waiting for SWISS-MODEL to build the loops (this can take 1-3 minutes)...")
        while not is_finished:
            status_response = requests.get(status_url, headers={"Authorization": f"Token {SM_TOKEN}"})
            status_data = status_response.json()
            status = status_data.get("status", "").upper()
            
            if status in ["COMPLETED", "SUCCESS", "DONE"]:
                is_finished = True
            elif status in ["FAILED", "ERROR"]:
                raise Exception("SWISS-MODEL processing engine failed to model this template.")
            else:
                print("Modeling execution in progress... checking again in 20 seconds.")
                time.sleep(20)
                
        # --- D. Download the Repaired PDB ---
        print("Modeling complete! Fetching repaired structural coordinates...")
        
        models = status_data.get("models", [])
        if not models:
            raise Exception("No generated structural coordinate models found.")
            
        download_url = models[0].get("coordinates_url")
        download_response = requests.get(download_url, headers={"Authorization": f"Token {SM_TOKEN}"})
        
        # Decompress the gzipped data from the server response
        try:
            decompressed_content = gzip.decompress(download_response.content)
        except Exception:
            # If the server returned it uncompressed, use it raw
            decompressed_content = download_response.content

        with open(repaired_pdb, "wb") as f:
            f.write(decompressed_content)
            
        print(f"[SUCCESS] Repaired receptor saved to: {repaired_pdb}")

        # =======================================================
        # NEW STEP: AUTOMATIC NATIVE RESIDUE NUMBERING RECOVERY
        # =======================================================
        print("\n--- RESTORING ORIGINAL NATIVE RESIDUE NUMBERING ---")
        
        # Reload both structures to map coordinates
        clean_struct = parser.get_structure("clean", cleaned_pdb)
        repaired_struct = parser.get_structure("repaired", repaired_pdb)
        
        # Gather original Alpha-Carbon locations and their true residue numbers
        orig_ca_mapping = []
        for model in clean_struct:
            for chain in model:
                for res in chain:
                    if 'CA' in res:
                        orig_ca_mapping.append((res['CA'].coord, chain.id, res.id[1]))
                        
        # Reposition and sync IDs on the repaired layout
        for model in repaired_struct:
            for chain in model:
                for res in chain:
                    if 'CA' in res:
                        rep_coord = res['CA'].coord
                        # Check distance to the closest original atom
                        distances = [np.linalg.norm(rep_coord - orig[0]) for orig in orig_ca_mapping]
                        if distances:
                            min_idx = np.argmin(distances)
                            # If they closely match spatially (< 0.5 Angstroms), map original ID
                            if distances[min_idx] < 0.5:
                                original_native_num = orig_ca_mapping[min_idx][2]
                                res.id = (res.id[0], original_native_num, res.id[2])

        # Overwrite file with correct native indices
        io_restore = PDBIO()
        io_restore.set_structure(repaired_struct)
        io_restore.save(repaired_pdb)
        print(f"[SUCCESS] Native PDB index formatting restored onto: {repaired_pdb}")

        # =======================================================
        # NEW STEP: RE-GRAFT STRIPPED WATERS & COFACTORS FROM CLEAN PDB
        # =======================================================
        print("\n--- RE-GRAFTING PRESERVED WATERS INTO REPAIRED RECEPTOR ---")
        waters_to_graft = []
        with open(cleaned_pdb, 'r') as cf:
            for line in cf:
                if line.startswith("HETATM") and "HOH" in line:
                    waters_to_graft.append(line)
                    
        if waters_to_graft:
            with open(repaired_pdb, 'r') as rf:
                repaired_lines = [line for line in rf if not line.startswith("END")]
            with open(repaired_pdb, 'w') as wf:
                wf.writelines(repaired_lines)
                wf.writelines(waters_to_graft)
                wf.write("END\n")
            print(f"[SUCCESS] Re-grafted {len(waters_to_graft)} water lines back into {repaired_pdb}")
        else:
            print("[INFO] No water molecules needed re-grafting.")

        final_receptor_name = repaired_pdb

    except Exception as e:
        print(f"\n[WARNING] SWISS-MODEL automated repair failed: {e}")
        print("Falling back to un-repaired clean PDB for config.")
        final_receptor_name = cleaned_pdb

    # =======================================================
    # NEW STEP: CONVERT TO PDBQT (via local MGLTools)
    # =======================================================
    print(f"\n--- CONVERTING {final_receptor_name} TO PDBQT ---")
    output_pdbqt = f"{final_receptor_name.split('.')[0]}.pdbqt"
    
    # Point directly to the double-nested MGLTools folder in your current directory
    pythonsh_path = "./mgltools_x86_64Linux2_1.5.7p1/mgltools_x86_64Linux2_1.5.7/bin/pythonsh"
    prep_script = "./mgltools_x86_64Linux2_1.5.7p1/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
    
    # Dynamically determine cleanup flags based on whether waters were grafted/kept
    if 'waters_to_graft' in locals() and len(waters_to_graft) > 0:
        cleanup_flags = "nphs_lps"        # Preserve kept water molecules
        print("[INFO] Preserving specific water molecules during PDBQT conversion.")
    else:
        cleanup_flags = "nphs_lps_waters" # Strip all waters if none were chosen
        print("[INFO] Cleaning up all water molecules during PDBQT conversion.")

    mgl_cmd = [
        pythonsh_path, 
        prep_script, 
        "-r", final_receptor_name,
        "-o", output_pdbqt,
        "-A", "checkhydrogens", 
        "-U", cleanup_flags
    ]

    try:
        print("Running local MGLTools to assign Kollman charges and add polar hydrogens...")
        
        # Run the command
        result = subprocess.run(mgl_cmd, capture_output=True, text=True, check=True)
        
        print(f"[SUCCESS] Receptor fully prepped and saved to: {output_pdbqt}")
        
    except FileNotFoundError:
        print(f"\n[ERROR] Could not find MGLTools executables at the specified paths.")
        print("Ensure the 'mgltools_x86_64Linux2_1.5.7p1' folder is fully extracted in the same directory as this script.")
        
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] MGLTools conversion failed. Error details:\n{e.stderr}\n{e.stdout}")

if __name__ == "__main__":
    main()

```
</details>
