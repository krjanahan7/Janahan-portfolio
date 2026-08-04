### 📌 Symmetry-Aware Ligand RMSD Calculation & PyMOL Superimposition

[![View Script](https://img.shields.io/badge/VIEW_SCRIPT-AUTO__RMSD__RENDER__ADVANCED.PY-blue?style=for-the-badge&logo=python)](#)
[![Status](https://img.shields.io/badge/STATUS-EXECUTABLE-brightgreen?style=for-the-badge)](#)

<details>
<summary>🔍 <b>Click to View Script Details & Usage Header</b></summary>

```python
# ==============================================================================
# SCRIPT NAME : auto_rmsd_render_advanced.py
# AUTHOR      : K. R. Janahan
# DESCRIPTION : Automatically scans PDB structures for ligands/ions, allows 
#               interactive selection, computes heavy-atom symmetry-aware RMSD 
#               via RDKit MCS, and generates a publication-quality PyMOL 3D 
#               superimposition rendering.
# DEPENDENCIES: PyMOL, RDKit, Python3
# USAGE       : python3 auto_rmsd_render_advanced.py
# ==============================================================================
```
</details>

<details>
<summary>💻 <b>Click to View / Copy Full Source Code</b></summary>

```python

import os
import sys
import pymol
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdFMCS

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

# CONFIGURATION
RAW_PDB = os.path.abspath("1_UNK-ORG.pdb")
COMPLEX_PDB = os.path.abspath("2_UNK-DOCK.pdb")
OUTPUT_IMAGE = "superimposed_result.png"

# Launch PyMOL headless (no GUI) early
pymol.finish_launching(["pymol", "-qc"])

STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HOH", "WAT", "SOL", "DOD"
}

def safe_load_pdb(pdb_path, obj_name):
    """
    Safely loads a PDB file into PyMOL. If Open Babel corrupted the CONECT
    records, this strips them out and loads clean ATOM/HETATM coordinates.
    """
    if not os.path.exists(pdb_path) or os.path.getsize(pdb_path) == 0:
        print(f"Error: File '{pdb_path}' not found or empty.")
        sys.exit(1)

    # First attempt standard loading
    cmd.load(pdb_path, obj_name)

    # If PyMOL rejected the file (common with Open Babel CONECT bugs), clean it:
    if obj_name not in cmd.get_object_list("all"):
        with open(pdb_path, "r", errors="ignore") as f:
            clean_lines = [
                line for line in f
                if line.startswith(("ATOM", "HETATM", "MODEL", "ENDMDL", "CRYST1"))
            ]
        cmd.read_pdbstr("".join(clean_lines), obj_name)

    # Final verification
    if obj_name not in cmd.get_object_list("all"):
        print(f"Fatal Error: PyMOL could not parse coordinate data from '{pdb_path}'.")
        sys.exit(1)


def scan_ligands_and_ions(pdb_file, obj_name):
    """Scans for non-standard residues across both ATOM and HETATM records."""
    safe_load_pdb(pdb_file, obj_name)
    het_residues = {}

    def collect_res(resn):
        clean_res = resn.strip().upper()
        if clean_res not in STANDARD_AMINO_ACIDS:
            het_residues[clean_res] = het_residues.get(clean_res, 0) + 1

    cmd.iterate(
        obj_name,
        "collect_res(resn)",
        space={"collect_res": collect_res, "STANDARD_AMINO_ACIDS": STANDARD_AMINO_ACIDS}
    )
    return het_residues


# --- STEP 0: Interactive Dual Residue Selection ---
print("=" * 60)
print("Scanning PDB files for ligands, cofactors, and ions...")
print("=" * 60)

cmd.reinitialize()
raw_het = scan_ligands_and_ions(RAW_PDB, "raw_scan")
complex_het = scan_ligands_and_ions(COMPLEX_PDB, "complex_scan")

if not raw_het and not complex_het:
    print("No ligands, cofactors, or ions found in either file.")
    cmd.quit()
    sys.exit(1)

COMMON_IONS = {"CL", "NA", "MG", "ZN", "CA", "FE", "K", "MN", "CU", "BR", "I"}

# Sort only the residues present inside the RAW_PDB
raw_resnames = sorted(list(raw_het.keys()))

print(f"\nFound {len(raw_resnames)} non-solvent residue(s) inside your reference file ({os.path.basename(RAW_PDB)}):\n")
print(f"{'ResName':<10} {'Type Guess':<20} {'Raw Atoms':<12}")
print("-" * 45)

for resn in raw_resnames:
    res_type = "Ion" if resn.strip() in COMMON_IONS else "Ligand / Cofactor"
    raw_count = raw_het[resn]
    print(f"{resn:<10} {res_type:<20} {raw_count:<12}")

print("-" * 45)

# 1. Ask for Reference Ligand Name
while True:
    choice_ref = input("\n[1/2] Enter the ResName for the REFERENCE crystal ligand: ").strip().upper()
    if choice_ref in raw_het:
        REF_RESNAME = choice_ref
        break
    else:
        print(f"'{choice_ref}' was not found in the reference file. Please enter one of: {', '.join(raw_resnames)}")

# 2. Ask for Docked Ligand Name
print("\n" + "-" * 50)
print(f"Available non-protein residues detected inside {os.path.basename(COMPLEX_PDB)}:")
for resn, count in complex_het.items():
    print(f"  -> {resn} ({count} atoms)")
print("-" * 50)

while True:
    choice_dock = input(f"[2/2] Enter the ResName for the DOCKED ligand : ").strip().upper()
    if not choice_dock:
        choice_dock = REF_RESNAME  # Default to the reference name if user presses Enter

    if choice_dock in complex_het:
        DOCKED_RESNAME = choice_dock
        break
    else:
        avail_options = ", ".join(complex_het.keys())
        print(f"'{choice_dock}' was not found inside {os.path.basename(COMPLEX_PDB)}. Please choose from: {avail_options}")

# --- STEP 1: PyMOL Extraction ---
print(f"\nExtracting reference ligand '{REF_RESNAME}' and docked ligand '{DOCKED_RESNAME}' using PyMOL...")

cmd.reinitialize()
safe_load_pdb(RAW_PDB, "raw")
safe_load_pdb(COMPLEX_PDB, "complex")

cmd.select("ref_lig", f"raw and resn {REF_RESNAME}")
cmd.select("docked_lig", f"complex and resn {DOCKED_RESNAME}")

cmd.save("ref_ligand.pdb", "ref_lig")
cmd.save("docked_ligand.pdb", "docked_lig")
cmd.reinitialize()

# --- STEP 2: RDKit Alignment & RMSD ---
print("Calculating symmetry-aware RMSD with RDKit...")

ref_mol = Chem.MolFromPDBFile("ref_ligand.pdb", removeHs=True)
docked_mol = Chem.MolFromPDBFile("docked_ligand.pdb", removeHs=True)

if not ref_mol or not docked_mol:
    print("Error: RDKit could not parse one of the extracted ligand files.")
    cmd.quit()
    sys.exit(1)

# Find common backbone/elements ignoring slight bond order guessing differences
mcs = rdFMCS.FindMCS(
    [ref_mol, docked_mol],
    atomCompare=rdFMCS.AtomCompare.CompareElements,
    bondCompare=rdFMCS.BondCompare.CompareAny,
    timeout=10,
)

patt = Chem.MolFromSmarts(mcs.smartsString)

# Get all symmetric atom mappings between reference and docked ligand
ref_matches = ref_mol.GetSubstructMatches(patt)
docked_matches = docked_mol.GetSubstructMatches(patt)

if not ref_matches or not docked_matches:
    print("Error: Could not establish atom mapping between ligands.")
    cmd.quit()
    sys.exit(1)

# Evaluate all symmetry permutations to find the optimal minimum RMSD
best_rmsd = float("inf")
best_alignment_map = None

for r_match in ref_matches:
    for d_match in docked_matches:
        atom_map = list(zip(d_match, r_match))
        try:
            val = rdMolAlign.AlignMol(docked_mol, ref_mol, atomMap=atom_map)
            if val < best_rmsd:
                best_rmsd = val
                best_alignment_map = atom_map
        except Exception:
            continue

rmsd = best_rmsd
# Perform final alignment using the best symmetry map
rdMolAlign.AlignMol(docked_mol, ref_mol, atomMap=best_alignment_map)

print(f"--> Heavy-Atom Symmetry-Aware RMSD: {rmsd:.3f} Å")

Chem.MolToPDBFile(docked_mol, "docked_ligand_aligned.pdb")

# --- STEP 3: PyMOL Rendering & Export ---
print("Rendering final high-resolution image...")
cmd.reinitialize()

# Load cleaned coordinate structures into the scene
safe_load_pdb("ref_ligand.pdb", "Reference")
safe_load_pdb("docked_ligand_aligned.pdb", "Docked")

# Save BOTH aligned ligands into one combined PDB file
COMBINED_PDB = "superimposed_pair.pdb"
cmd.save(COMBINED_PDB, "Reference or Docked")
print(f"--> Saved combined 3D coordinates to '{COMBINED_PDB}'")

# Reset representations
cmd.hide("everything")
cmd.show("sticks", "Reference or Docked")

# Apply clean carbon-only coloring so heteroatoms (O/N/S) stay distinct
cmd.util.cbag("Reference")  # Cyan carbons for reference
cmd.util.cbam("Docked")     # Magenta carbons for docked

# Center view on the superimposed pair
cmd.center("Reference")
cmd.orient("Reference")

# Increase buffer significantly (6.5) to give plenty of white canvas space below the molecule
cmd.zoom("Reference or Docked", buffer=6.5)

# Create anchor directly on the reference molecule
cmd.pseudoatom("rmsd_lbl", selection="Reference")
cmd.label("rmsd_lbl", f'"LRMSD = {rmsd:.2f} A"')
cmd.show("label", "rmsd_lbl")

# Format label text properties
cmd.set("label_color", "black", "rmsd_lbl")
cmd.set("label_size", 28, "rmsd_lbl")
cmd.set("label_font_id", 7, "rmsd_lbl")

# Set Y to -22.0 to push the label way down below the structure
cmd.set("label_position", [0.0, -22.0, 2.0], "rmsd_lbl")

# Visual polish for rendering
cmd.set("stick_radius", 0.22)
cmd.set("ray_opaque_background", 1)
cmd.bg_color("white")

# Ray trace and export
cmd.ray(2000, 2000)
cmd.png(OUTPUT_IMAGE, dpi=300)

print(f"\nSuccess! Saved clean rendering to '{OUTPUT_IMAGE}'")
cmd.quit()

```
</details>
