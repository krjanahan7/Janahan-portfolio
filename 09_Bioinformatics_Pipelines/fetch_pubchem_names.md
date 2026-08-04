### 📌 PubChem Chemical Identifier to Compound Name Converter[cite: 6]

[![View Script](https://img.shields.io/badge/VIEW_SCRIPT-FETCH__PUBCHEM__NAMES.PY-blue?style=for-the-badge&logo=python)](#)
[![Status](https://img.shields.io/badge/STATUS-EXECUTABLE-brightgreen?style=for-the-badge)](#)

<details>
<summary>🔍 <b>Click to View Script Details & Usage Header</b></summary>

```python
# ==============================================================================
# SCRIPT NAME : fetch_pubchem_names.py
# AUTHOR      : K. R. Janahan
# DESCRIPTION : Fetches IUPAC/common compound names from PubChem PUG-REST API 
#               using CIDs and appends them to a CSV screening report. Includes 
#               rate-limiting pauses to ensure API compliance.
# DEPENDENCIES: Python3 (Standard Libraries: csv, urllib, time)
# USAGE       : python3 fetch_pubchem_names.py

# ==============================================================================
```
</details>

<details>
<summary>💻 <b>Click to View / Copy Full Source Code</b></summary>

```python
import csv
import urllib.request
import time
import sys
import os

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

input_file = 'Final_Screening_Report.csv'
output_file = 'Final_Screening_Report_with_Names.csv'

def get_compound_name(cid):
    # PubChem's API endpoint that returns just the Title as plain text
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Title/TXT"
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            return response.read().decode('utf-8').strip()
    except urllib.error.URLError:
        return "Name_Not_Found"

print("Starting PubChem API lookup...")

try:
    with open(input_file, mode='r') as infile, open(output_file, mode='w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        # Read the header and add the new column
        header = next(reader)
        header.append("Compound_Name")
        writer.writerow(header)
        
        # Identify which column contains the CID
        try:
            cid_index = header.index("CID")
        except ValueError:
            print("Error: Could not find a column named 'CID' in your CSV header.")
            sys.exit(1)
            
        # Process each row
        count = 0
        for row in reader:
            if not row: continue  # Skip empty lines
            cid = row[cid_index]
            
            print(f"Fetching name for CID {cid}...", end=" ")
            name = get_compound_name(cid)
            print(f"Found: {name}")
            
            row.append(name)
            writer.writerow(row)
            count += 1
            
            # PubChem restricts traffic to 5 requests per second. 
            # This tiny pause prevents your IP from getting temporarily banned.
            time.sleep(0.25)

    print(f"\n✅ Success! Processed {count} molecules.")
    print(f"📂 Updated report saved as: {output_file}")

except FileNotFoundError:
    print(f"Error: Could not find {input_file} in this directory.")

```
</details>
