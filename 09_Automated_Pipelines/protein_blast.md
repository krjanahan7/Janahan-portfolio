### 📌 PDB Chain Sequence Extraction & RCSB Structural Similarity Search[cite: 11]

[![View Script](https://img.shields.io/badge/VIEW_SCRIPT-PROTEIN_BLAST.PY-blue?style=for-the-badge&logo=python)](#)
[![Status](https://img.shields.io/badge/STATUS-EXECUTABLE-brightgreen?style=for-the-badge)](#)

<details>
<summary>🔍 <b>Click to View Script Details & Usage Header</b></summary>

```python
# ==============================================================================
# SCRIPT NAME : protein_blast.py
# AUTHOR      : K. R. Janahan
# DESCRIPTION : Extracts FASTA sequences from PDB chain coordinates, queries the 
#               RCSB PDB Sequence Similarity Search API (v2), and generates a 
#               ranked report of matching structural PDB IDs with identity 
#               percentages and E-values.
# DEPENDENCIES: Biopython, Requests, Python3
# USAGE       : python3 protein_blast.py
# ==============================================================================
```
</details>

<details>
<summary>💻 <b>Click to View / Copy Full Source Code</b></summary>

```python

import sys
import subprocess

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

# --- AUTOMATIC DEPENDENCY CHECKER & INSTALLER ---
REQUIRED_PACKAGES = {
    "Bio": "biopython",
    "requests": "requests"
}

def install_and_import_dependencies():
    """Checks for required packages and installs them if missing."""
    missing_packages = []
    for module_name, pip_name in REQUIRED_PACKAGES.items():
        try:
            __import__(module_name)
        except ImportError:
            missing_packages.append(pip_name)
            
    if missing_packages:
        print(f"Missing required packages: {', '.join(missing_packages)}")
        print("Attempting to install them automatically...")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", *missing_packages])
            print("Successfully installed missing dependencies.\n")
        except subprocess.CalledProcessError as e:
            print(f"\nError: Automated installation failed with exit code {e.returncode}.")
            print("Please try installing the packages manually by running:")
            print(f"pip install {' '.join(missing_packages)}")
            sys.exit(1)

# Run the installation checker
install_and_import_dependencies()

# --- MAIN SCRIPT IMPORTS ---
import os
import warnings
from Bio import SeqIO
import requests

# Suppress Biopython parser & construction warnings to keep terminal clean
warnings.filterwarnings("ignore")


# --- PIPELINE FUNCTIONS ---
def extract_fasta_from_pdb(pdb_filepath):
    """Extracts FASTA sequences for all protein chains in a PDB file."""
    sequences = []
    try:
        for record in SeqIO.parse(pdb_filepath, "pdb-atom"):
            sequences.append((record.id, str(record.seq)))
    except Exception as e:
        print(f"Error parsing PDB file: {e}")
    return sequences


def search_rcsb_pdb(sequence, min_identity=0.4):
    """
    Queries RCSB PDB Sequence Similarity Search API (v2) for matching PDB IDs,
    extracting alignment identity percentages and E-values.
    """
    print("Querying RCSB PDB for top matching structural IDs...")
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    query_payload = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 0.1,
                "identity_cutoff": min_identity,
                "sequence_type": "protein",
                "value": sequence
            }
        },
        "request_options": {
            "scoring_strategy": "sequence"
        },
        "return_type": "polymer_entity"
    }
    
    try:
        response = requests.post(url, json=query_payload, timeout=30)
        if response.status_code == 200:
            results = response.json().get("result_set", [])
            pdb_matches = []
            seen_ids = set()
            
            for item in results:
                identifier = item.get("identifier", "")
                pdb_id = identifier.split("_")[0]
                
                # Ensure we only keep unique 4-character PDB IDs
                if pdb_id in seen_ids:
                    continue
                seen_ids.add(pdb_id)
                
                # Safely extract alignment statistics from the API response
                evalue_str = "N/A"
                identity_str = "N/A"
                try:
                    match_context = item["services"][0]["nodes"][0]["match_context"][0]
                    
                    evalue_val = match_context.get("evalue")
                    if isinstance(evalue_val, (int, float)):
                        evalue_str = f"{evalue_val:.2e}"
                    else:
                        evalue_str = str(evalue_val)
                        
                    ident_val = match_context.get("sequence_identity")
                    if isinstance(ident_val, (int, float)):
                        identity_str = f"{ident_val * 100:.1f}%"
                    else:
                        identity_str = str(ident_val)
                except (KeyError, IndexError, TypeError):
                    # Fallback if specific alignment context is unavailable
                    score_val = item.get("score", "N/A")
                    identity_str = f"Score: {score_val}"
                    
                pdb_matches.append({
                    "pdb_id": pdb_id,
                    "identity": identity_str,
                    "evalue": evalue_str
                })
                
                if len(pdb_matches) >= 10:
                    break
                    
            return pdb_matches
        else:
            print(f"RCSB Search API returned status code: {response.status_code}")
            return []
    except Exception as e:
        print(f"Error connecting to RCSB API: {e}")
        return []


def main(pdb_file):
    if not os.path.exists(pdb_file):
        print(f"Error: File '{pdb_file}' not found.")
        return

    chains = extract_fasta_from_pdb(pdb_file)
    if not chains:
        print("No valid polypeptide chains found in the PDB file.")
        return

    for chain_id, seq in chains:
        print(f"\n=========================================================================")
        print(f"Processing Chain: {chain_id} (Length: {len(seq)} aa)")
        print(f"Sequence preview: {seq[:45]}...")
        print(f"=========================================================================")
        
        # Pull Top PDB matches directly from RCSB PDB
        top_matches = search_rcsb_pdb(seq)
        
        print("\n--- Top Matching PDB IDs (RCSB PDB Website) ---")
        if top_matches:
            # Print table header
            print(f"{'Rank':<6}{'PDB ID':<8}{'Identity':<12}{'E-value':<12}{'Structure Link'}")
            print("-" * 75)
            
            # Print ranked rows
            for idx, match in enumerate(top_matches, 1):
                pid = match['pdb_id']
                ident = match['identity']
                eval_str = match['evalue']
                link = f"https://www.rcsb.org/structure/{pid}"
                print(f"{idx:<6}{pid:<8}{ident:<12}{eval_str:<12}{link}")
            print("-" * 75)
        else:
            print("No significant PDB structure matches found.")


if __name__ == "__main__":
    # Dynamically check if you provided a file in the terminal command
    if len(sys.argv) > 1:
        PDB_FILE_PATH = sys.argv[1]
    else:
        PDB_FILE_PATH = "unk.pdb"
        print(f"No file specified. Defaulting to '{PDB_FILE_PATH}'")
        
    main(PDB_FILE_PATH)
```
</details>
