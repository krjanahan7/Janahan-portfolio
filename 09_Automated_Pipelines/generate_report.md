### 📌 AutoDock Vina Logs to PubChem Multi-Property Report Generator[cite: 7]

[![View Script](https://img.shields.io/badge/VIEW_SCRIPT-EXTRACT__VINA__PROPERTIES.PY-blue?style=for-the-badge&logo=python)](#)
[![Status](https://img.shields.io/badge/STATUS-EXECUTABLE-brightgreen?style=for-the-badge)](#)

<details>
<summary>🔍 <b>Click to View Script Details & Usage Header</b></summary>

```python
# ==============================================================================
# SCRIPT NAME : generate_report.py
# AUTHOR      : K. R. Janahan
# DESCRIPTION : Direct pipeline that extracts AutoDock Vina top-pose binding 
#               affinities from log files, queries PubChem API for molecular 
#               properties and IUPAC titles, and writes a unified CSV report.
# DEPENDENCIES: Python3 (Standard Libraries: os, csv, glob, urllib, json, time)
# USAGE       : python3 generate_report.py

# ==============================================================================
```
</details>

<details>
<summary>💻 <b>Click to View / Copy Full Source Code</b></summary>

```python

import os
import csv
import glob
import urllib.request
import json
import time
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

# --- Configuration ---
final_output_file = 'Screening_Report.csv'

def get_compound_data(cid):
    # Fetch Title (Name), MolecularWeight, XLogP, and TPSA in a single JSON request
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Title,MolecularWeight,XLogP,TPSA/JSON"
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            raw_data = json.loads(response.read().decode('utf-8'))
            properties = raw_data['PropertyTable']['Properties'][0]
            
            # Use .get() to gracefully handle missing data
            name = properties.get('Title', 'Name_Not_Found')
            mol_wt = properties.get('MolecularWeight', 'N/A')
            logp = properties.get('XLogP', 'N/A')
            tpsa = properties.get('TPSA', 'N/A')
            
            return name, mol_wt, logp, tpsa
    except Exception:
        return "Error", "Error", "Error", "Error"

print("Starting direct pipeline: Vina Logs -> PubChem -> Final CSV")

log_files = glob.glob('CID_*/log.txt')
if not log_files:
    print("Error: No 'CID_*/log.txt' files found in the current directory.")
    exit()

# Open the final CSV file for writing row-by-row
with open(final_output_file, mode='w', newline='') as outfile:
    writer = csv.writer(outfile)
    
    # Write the comprehensive header 
    writer.writerow(["CID", "BindingAffinity(kcal/mol)", "MolecularWeight", "LogP", "TPSA", "Compound_Name"])
    
    count = 0
    for log_path in log_files:
        # 1. Extract the CID from the folder name
        folder_name = os.path.dirname(log_path)
        cid = folder_name.replace('CID_', '')
        
        # 2. Extract Binding Affinity
        binding_affinity = "N/A"
        with open(log_path, 'r') as file:
            for line in file:
                parts = line.split()
                # The top pose in AutoDock Vina always starts with '1' under the mode column
                if len(parts) > 1 and parts[0] == '1':
                    binding_affinity = parts[1]
                    break 
                    
        # 3. Fetch PubChem Data
        print(f"Processing CID {cid}...", end=" ")
        name, mol_wt, logp, tpsa = get_compound_data(cid)
        
        # Print a concise summary to the console
        display_name = name[:20] + '...' if len(name) > 20 else name
        print(f"Affinity: {binding_affinity} | MW: {mol_wt} | LogP: {logp} | [{display_name}]")
        
        # 4. Write directly to the final CSV
        writer.writerow([cid, binding_affinity, mol_wt, logp, tpsa, name])
        count += 1
        
        # 5. Rest for 0.25s to respect PubChem's rate limit
        time.sleep(0.25)

print(f"\n✅ Success! Processed {count} molecules in a single pass.")
print(f"📂 Final report saved as: {final_output_file}")

```
</details>
