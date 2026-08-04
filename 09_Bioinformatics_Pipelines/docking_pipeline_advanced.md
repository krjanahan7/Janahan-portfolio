### 📌 High-Throughput Virtual Screening & Docking Pipeline (AutoDock Vina / OpenBabel)[cite: 5]

[![View Script](https://img.shields.io/badge/VIEW_SCRIPT-VIRTUAL__SCREENING__PIPELINE.PY-blue?style=for-the-badge&logo=python)](#)
[![Status](https://img.shields.io/badge/STATUS-EXECUTABLE-brightgreen?style=for-the-badge)](#)

<details>
<summary>🔍 <b>Click to View Script Details & Usage Header</b></summary>

```python
# ==============================================================================
# SCRIPT NAME : docking_pipeline_advanced.py
# AUTHOR      : K. R. Janahan
# DESCRIPTION : End-to-end interactive virtual screening pipeline. Performs 
#               Lipinski filtering, OpenBabel format conversion, AutoDock Vina 
#               docking simulations, automated report generation, and protein-ligand 
#               complex structure construction.
# DEPENDENCIES: AutoDock Vina, OpenBabel, Pandas, Requests, Python3
# USAGE       : python3 docking_pipeline_advanced.py

# ==============================================================================
```
</details>

<details>
<summary>💻 <b>Click to View / Copy Full Source Code</b></summary>

```python

import os
import time
import subprocess
import requests
import pandas as pd
import re
import shutil
import sys

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

# ==========================================
# 1. SETTINGS & PARAMETERS
# ==========================================
INPUT_FOLDER = "Penicillin_G_90_similarity"  
MAX_MOLECULES = 1                            

RECEPTOR = "5AXO.pdbqt"
CONF_FILE = "./Autodock-Vina/conf.txt"
PUG_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

for folder in ["pdbqt", "results", "reports"]:
    os.makedirs(folder, exist_ok=True)

session = requests.Session()
session.headers.update({"User-Agent": "Bioinformatics-Screening-Script/3.0"})

compound_properties = {}
valid_cids = []
valid_files = [] 

# ==========================================
# 2. INTERACTIVE MENU
# ==========================================
print("\n" + "="*50)
print(" 🧬 VIRTUAL SCREENING PIPELINE")
print("="*50)
print("1. Scan Local Folder & Apply Lipinski Filter")
print("2. Convert SDF to PDBQT (OpenBabel)")
print("3. Molecular Docking Simulations (Vina)")
print("4. Organize Results Directory")
print("5. Metrics Compilation, Name Fetching & Reporting")
print("6. Extract Poses, Convert to PDB & Merge")
print("="*50)

user_input = input("\nEnter the steps you want to run (e.g., 1,2,3 or 4,5,6) or 'all': ").strip().lower()

if user_input == 'all':
    run_steps = [1, 2, 3, 4, 5, 6]
else:
    # Parse the comma-separated numbers into a list of integers
    run_steps = [int(s.strip()) for s in user_input.split(",") if s.strip().isdigit()]

if not run_steps:
    print("⚠️ No valid steps selected. Exiting.")
    exit()

# ==========================================
# 3. FILTER FUNCTIONS
# ==========================================
def fetch_and_verify_properties(cid):
    try:
        url = (
            f"{PUG_URL}/compound/cid/{cid}/property/"
            "MolecularWeight,XLogP,HBondDonorCount,"
            "HBondAcceptorCount,TPSA/JSON"
        )
        response = session.get(url, timeout=30).json()
        props = response["PropertyTable"]["Properties"][0]

        mw = float(props.get("MolecularWeight", 0))
        logp = float(props.get("XLogP", 0))
        hbd = int(props.get("HBondDonorCount", 0))
        hba = int(props.get("HBondAcceptorCount", 0))
        tpsa = float(props.get("TPSA", 0)) 

        compound_properties[str(cid)] = {
            "MW": mw, "LogP": logp, "HBD": hbd, "HBA": hba, "TPSA": tpsa
        }

        lipinski_pass = (mw <= 500) and (logp <= 5) and (hbd <= 5) and (hba <= 10)
        return lipinski_pass
    except Exception:
        return True 

# ==========================================
# PIPELINE EXECUTION STEPS
# ==========================================

# STEP 1
if 1 in run_steps:
    print(f"\n[STEP 1] Scanning local folder '{INPUT_FOLDER}' & Applying Lipinski filter...")
    if not os.path.exists(INPUT_FOLDER):
        print(f"⚠️ Critical Error: The folder '{INPUT_FOLDER}' was not found.")
        exit(1)

    sdf_files = sorted([f for f in os.listdir(INPUT_FOLDER) if f.endswith(".sdf")])
    print(f"-> Found {len(sdf_files)} SDF files. Limit set to {MAX_MOLECULES}.")

    for idx, filename in enumerate(sdf_files, 1):
        if len(valid_cids) >= MAX_MOLECULES:
            print(f"\n   🛑 Reached limit of {MAX_MOLECULES} molecules.")
            break

        cid_match = re.search(r'\d+', filename)
        cid_candidate = cid_match.group() if cid_match else f"unknown_{idx}"

        print(f" Processing [{idx}/{len(sdf_files)}] | File: {filename} | CID: {cid_candidate}")
        
        if str(cid_candidate).isdigit():
            if not fetch_and_verify_properties(cid_candidate):
                print("   ✖ Excluded: Fails Lipinski's Rule of 5.")
                continue
            time.sleep(0.2) 
        else:
            print("   ⚠️ No CID found. Proceeding without Lipinski filter.")

        valid_cids.append(cid_candidate)
        valid_files.append((filename, cid_candidate))

# STEP 2
if 2 in run_steps:
    print("\n[STEP 2] Converting structural formats (SDF to PDBQT via OpenBabel)...")
    
    # Safeguard: If Step 1 was skipped, populate valid_files from the folder manually
    if not valid_files:
        sdf_files = sorted([f for f in os.listdir(INPUT_FOLDER) if f.endswith(".sdf")])
        for idx, filename in enumerate(sdf_files, 1):
            cid_match = re.search(r'\d+', filename)
            cid_candidate = cid_match.group() if cid_match else f"unknown_{idx}"
            valid_files.append((filename, cid_candidate))

    for idx, (filename, cid_extracted) in enumerate(valid_files, 1):
        input_sdf = os.path.join(INPUT_FOLDER, filename)
        output_pdbqt = os.path.join("pdbqt", f"FINAL-PDBQT_{cid_extracted}.pdbqt")

        if os.path.exists(output_pdbqt):
            print(f" [{idx}/{len(valid_files)}] {filename} -> Output already exists. Skipping.")
            continue

        print(f" [{idx}/{len(valid_files)}] Converting: {filename} -> PDBQT")
        subprocess.run(["obabel", input_sdf, "-O", output_pdbqt, "-p", "7.4"], 
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

# STEP 3
if 3 in run_steps:
    print(f"\n[STEP 3] Initiating Molecular Docking against Target: {RECEPTOR}...")
    pdbqt_inventory = [f for f in os.listdir("pdbqt") if f.startswith("FINAL-PDBQT_")]

    if not os.path.exists(CONF_FILE):
        with open(CONF_FILE, "w") as f:
            f.write("# Auto-generated mock conf template\nsize_x = 20\nsize_y = 20\nsize_z = 20\n")

    for idx, ligand_file in enumerate(pdbqt_inventory, 1):
        ligand_path = os.path.join("pdbqt", ligand_file)
        ligand_id = ligand_file.replace("FINAL-PDBQT_", "").replace(".pdbqt", "")
        
        output_mesh = f"results/FINAL-PDBQT_{ligand_id}_out.pdbqt"
        log_dump = f"results/FINAL-PDBQT_{ligand_id}_log.txt"

        if os.path.exists(output_mesh):
            print(f" [{idx}/{len(pdbqt_inventory)}] Compound {ligand_id} already docked. Skipping.")
            continue

        print(f" [{idx}/{len(pdbqt_inventory)}] Running AutoDock Vina simulation for CID: {ligand_id}...")
        execution = subprocess.run([
            "./Autodock-Vina/vina", "--config", CONF_FILE, "--receptor", RECEPTOR,
            "--ligand", ligand_path, "--out", output_mesh
        ], capture_output=True, text=True)

        with open(log_dump, "w") as log:
            log.write(execution.stdout)
            log.write(execution.stderr)

# STEP 4 (Organize Results Directory)
if 4 in run_steps:
    print("\n[STEP 4] Organizing output files in 'results' directory...")
    logs_dir = os.path.join("results", "logs")
    pdbqt_dir = os.path.join("results", "out.pdbqt")
    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(pdbqt_dir, exist_ok=True)

    organized_count = 0
    if os.path.exists("results"):
        for file in os.listdir("results"):
            file_path = os.path.join("results", file)
            if os.path.isdir(file_path): continue
                
            if file.endswith("_log.txt"):
                shutil.move(file_path, os.path.join(logs_dir, file))
                organized_count += 1
            elif file.endswith("_out.pdbqt"):
                shutil.move(file_path, os.path.join(pdbqt_dir, file))
                organized_count += 1

    print(f" -> Successfully organized {organized_count} files.")

# STEP 5 (Merged: Metrics Compilation + Name Fetching + Reporting)
if 5 in run_steps:
    print("\n[STEP 5] Parsing results, fetching compound names, and compiling report...")
    import urllib.request
    
    def get_compound_name(cid):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Title/TXT"
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req) as response:
                return response.read().decode('utf-8').strip()
        except Exception:
            return "Name_Not_Found"

    dataset = []
    pdbqt_dir = os.path.join("results", "out.pdbqt")
    
    if os.path.exists(pdbqt_dir) and os.listdir(pdbqt_dir):
        # 1. Parse Docking Results
        for file in os.listdir(pdbqt_dir):
            if not file.endswith("_out.pdbqt"):
                continue

            filepath = os.path.join(pdbqt_dir, file)
            binding_energy = None

            with open(filepath, "r") as f:
                for line in f:
                    match = re.search(r"REMARK VINA RESULT:\s*(-?\d+\.\d+)", line)
                    if match:
                        binding_energy = float(match.group(1))
                        break

            if binding_energy is not None:
                cid_key = file.replace("FINAL-PDBQT_", "").replace("_out.pdbqt", "")
                meta = compound_properties.get(cid_key, {"MW": "N/A", "LogP": "N/A", "TPSA": "N/A"})
                dataset.append({
                    "CID": cid_key, 
                    "Compound_Name": "Pending...", # To be fetched
                    "BindingAffinity(kcal/mol)": binding_energy,
                    "MolecularWeight": meta["MW"], 
                    "LogP": meta["LogP"], 
                    "TPSA": meta["TPSA"]
                })

        # 2. Fetch Compound Names
        print(f" -> Fetching names for {len(dataset)} compounds from PubChem...")
        for row in dataset:
            print(f"    CID {row['CID']}...", end=" ", flush=True)
            name = get_compound_name(row["CID"])
            print(f"Found: {name}")
            row["Compound_Name"] = name
            time.sleep(0.25)

        # 3. Build DataFrame and Generate Reports
        df = pd.DataFrame(dataset)
        if not df.empty:
            # Reorder columns for readability
            df = df[["CID", "Compound_Name", "BindingAffinity(kcal/mol)", "MolecularWeight", "LogP", "TPSA"]]
            df = df.sort_values(by=["BindingAffinity(kcal/mol)", "MolecularWeight"], ascending=[True, True])
            
            report_csv = "reports/Full_Screening_Report.csv"
            df.to_csv(report_csv, index=False)
            
            best_compound, worst_compound = df.iloc[0], df.iloc[-1]
            good_targets = df[df["BindingAffinity(kcal/mol)"] <= -7.0]
            
            report_path = "reports/Executive_Summary.md"
            with open(report_path, "w") as rep:
                rep.write(f"# Virtual Screening Pipeline Report: Target {RECEPTOR}\n")
                rep.write(f"**Input Directory:** {INPUT_FOLDER}\n")
                rep.write(f"**Molecules Processed:** {len(dataset)}\n\n")
                rep.write("## 🏆 Lead Summary Metrics\n")
                rep.write(f"* **Best Compound:** {best_compound['Compound_Name']} (CID {best_compound['CID']}) -> {best_compound['BindingAffinity(kcal/mol)']} kcal/mol\n")
                rep.write(f"* **Worst Compound:** {worst_compound['Compound_Name']} (CID {worst_compound['CID']}) -> {worst_compound['BindingAffinity(kcal/mol)']} kcal/mol\n")
                rep.write(f"* **Promising Hits (<= -7.0 kcal/mol):** {len(good_targets)}\n\n")
                rep.write("## 📊 Top 10 Screening Matrix\n")
                rep.write(df.head(10).to_markdown(index=False))
                
            print("\n ✅ PROCESSING COMPLETE")
            print(f"    Saved Dataset: {report_csv}")
            print(f"    Saved Summary: {report_path}")
        else:
            print("❌ No valid binding criteria data could be aggregated.")
    else:
        print("⚠️ No docking results found in 'results/out.pdbqt' to report.")


# STEP 6
if 6 in run_steps:
    print("\n[STEP 6] Extracting Top Poses, Converting to PDB, and Generating Complexes...")
    logs_dir = os.path.join("results", "logs")
    pdbqt_dir = os.path.join("results", "out.pdbqt")

    def parse_vina_log_adapted(log_path):
        try:
            with open(log_path, 'r') as f:
                for line in f:
                    if line.strip().startswith("-----+"):
                        parts = next(f).split()
                        if len(parts) >= 3:
                            return float(parts[1]), float(parts[2])
        except Exception:
            pass
        return None, None

    def extract_first_pose(pdbqt_path, output_path):
        with open(pdbqt_path, 'r') as infile, open(output_path, 'w') as outfile:
            in_model = False
            for line in infile:
                if line.startswith("MODEL"): in_model = True
                if in_model: outfile.write(line)
                if line.startswith("ENDMDL") and in_model: break

    extract_output_dir = os.path.join("results", "complex-formation")
    lig_pdb_dir = os.path.join(extract_output_dir, "drugs.pdb")
    complex_dir = os.path.join(extract_output_dir, "complex.pdb")
    os.makedirs(lig_pdb_dir, exist_ok=True)
    os.makedirs(complex_dir, exist_ok=True)

    receptor_pdb = RECEPTOR.replace(".pdbqt", ".pdb")
    if not os.path.exists(receptor_pdb):
        subprocess.run(["obabel", "-ipdbqt", RECEPTOR, "-O", receptor_pdb, "-d"], stdout=subprocess.DEVNULL)

    extract_results = []
    seen_ligands = set()

    if os.path.exists(pdbqt_dir):
        for out_file in os.listdir(pdbqt_dir):
            if out_file.endswith("_out.pdbqt"):
                ligand_id = out_file.replace("FINAL-PDBQT_", "").replace("_out.pdbqt", "")
                if ligand_id in seen_ligands: continue
                seen_ligands.add(ligand_id)

                pdbqt_path = os.path.join(pdbqt_dir, out_file)
                log_path = os.path.join(logs_dir, f"FINAL-PDBQT_{ligand_id}_log.txt")
                
                docking_score, rmsd = parse_vina_log_adapted(log_path)
                if docking_score is not None:
                    extract_results.append([ligand_id, docking_score, rmsd])

                    first_pose = os.path.join(pdbqt_dir, f"{ligand_id}_first.pdbqt")
                    ligand_pdb = os.path.join(lig_pdb_dir, f"{ligand_id}_drug.pdb")
                    complex_pdb = os.path.join(complex_dir, f"{ligand_id}_complex.pdb")

                    extract_first_pose(pdbqt_path, first_pose)
                    subprocess.run(["obabel", "-ipdbqt", first_pose, "-O", ligand_pdb, "-d"], stdout=subprocess.DEVNULL)
                    
                    with open(receptor_pdb) as r, open(ligand_pdb) as l, open(complex_pdb, 'w') as out:
                        out.write(r.read() + l.read())

                    if os.path.exists(first_pose): os.remove(first_pose)

    print(f"\n🎉 Processed {len(extract_results)} hits successfully!")

```
</details>
