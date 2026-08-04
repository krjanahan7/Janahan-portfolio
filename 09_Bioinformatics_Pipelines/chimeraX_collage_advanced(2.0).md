### 📌 ChimeraX Batch H-Bond Analysis, Multi-View Rendering & Word Report Generator[cite: 10]

[![View Script](https://img.shields.io/badge/VIEW_SCRIPT-CHIMERAX__HBOND__BATCH__RENDER.PY-blue?style=for-the-badge&logo=python)](#)
[![Status](https://img.shields.io/badge/STATUS-EXECUTABLE-brightgreen?style=for-the-badge)](#)

<details>
<summary>🔍 <b>Click to View Script Details & Usage Header</b></summary>

```python
# ==============================================================================
# SCRIPT NAME : chimeraX_collage_advanced(2.0).py
# AUTHOR      : K. R. Janahan
# DESCRIPTION : ChimeraX automated batch processing script that calculates protein-ligand 
#               hydrogen bonds using dynamic geometric tolerances, exports full & zoomed 
#               high-resolution renders to sorted folders, and builds a unified Word (.docx) 
#               interaction report.
# DEPENDENCIES: ChimeraX, python-docx, Python3
# USAGE       : chimerax --nogui --offscreen chimeraX_collage_advanced(2.0).py
# ==============================================================================
```
</details>

<details>
<summary>💻 <b>Click to View / Copy Full Source Code</b></summary>

```python

# --- CHIMERAX 'ASK-ONCE' BATCH SCRIPT (H-BONDS ONLY, DOCX EXPORT, FOLDER SORTING & DYNAMIC SLOP) ---

from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure
import os
import glob 
import math 
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
# CONFIGURATION SETTINGS
# ==========================================
SEARCH_PATTERN = "*_complex.pdb"
# ==========================================

def process_batch(session):
    print("\n" + "="*50)
    print(" 1. GLOBAL BATCH SETUP")
    print("="*50)
    
    pdb_files = sorted(glob.glob(SEARCH_PATTERN))
    
    if not pdb_files:
        print(f"[ERROR] Could not find any files matching '{SEARCH_PATTERN}'.")
        return

    print(f"Found {len(pdb_files)} files. Inspecting the first file ({pdb_files[0]}) to gather setup options...\n")

    # --- SETUP SEPARATE OUTPUT DIRECTORIES ---
    full_view_dir = "Full_Renders"
    zoomed_view_dir = "Zoomed_Renders"
    os.makedirs(full_view_dir, exist_ok=True)
    os.makedirs(zoomed_view_dir, exist_ok=True)
    print(f"[INFO] Image destinations initialized: '{full_view_dir}/' and '{zoomed_view_dir}/'")

    std_amino_acids = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", 
                       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}

    # --- PHASE 1: INSPECT FIRST FILE & ASK QUESTIONS ONCE ---
    run(session, "close session")
    run(session, f"open {pdb_files[0]}")
    
    structures = session.models.list(type=AtomicStructure)
    if not structures:
        print("[ERROR] First file failed to load. Aborting.")
        return

    # Inspect Chains
    print("\n--- DETECTED PROTEIN CHAINS ---")
    chain_ids = set()
    for s in structures:
        for c in s.chains:
            if len([r for r in c.residues if r is not None]) > 0:
                cid = c.chain_id.strip() if c.chain_id and c.chain_id.strip() else "NO-ID"
                chain_ids.add(cid)
    
    print(f"Chains present in file: {', '.join(sorted(chain_ids))}")
    global_chains_input = input("Enter the chain ID(s) to keep for ALL files (comma-separated, or press Enter for ALL): ").strip().upper()

    # Inspect Ligands
    print("\n--- DETECTED LIGANDS/IONS ---")
    all_ligands = []
    for s in structures:
        for r in s.residues:
            if r is not None and r.name not in std_amino_acids and r.name != "HOH":
                all_ligands.append(r)
                
    if not all_ligands:
        print("None found.")
    else:
        for i, lig in enumerate(all_ligands):
            cid = lig.chain_id.strip() if lig.chain_id and lig.chain_id.strip() else "NO-ID"
            print(f"[{i}] {lig.name} (Model #{lig.structure.id_string}, Chain {cid}, Residue {lig.number})")
            
    global_ligands_input = input("\nEnter the numbers [ ] of the ligands to keep for ALL files (comma-separated, or press Enter for none): ").strip()

    # Render & Criteria Settings
    print("\n--- RENDER SETTINGS ---")
    global_zoom = input("Zoom factor (e.g., 0.85) [Press Enter for 0.85]: ").strip()
    if not global_zoom: global_zoom = "0.85"
    
    # --- UPGRADE: DYNAMIC HYDROGEN BONDING TOLERANCE PROMPT ---
    global_slop = input("Enter H-bond distance tolerance slop in Å (e.g., 0.4) [Press Enter for 0.4]: ").strip()
    if not global_slop: global_slop = "0.4"
    
    global_show_interact = input("Label protein-ligand Hydrogen Bonds for ALL files? (y/n): ").strip().lower()

    # --- AUTOMATIC DEPENDENCY CHECK & INITIALIZE WORD DOCX ---
    try:
        import docx
    except ImportError:
        print("[INFO] 'python-docx' library not found. Attempting automatic installation inside ChimeraX...")
        try:
            run(session, "pip install python-docx")
            import docx
        except Exception as e:
            print(f"[ERROR] Failed to install python-docx automatically: {e}")
            print("Please run the command 'pip install python-docx' manually in the ChimeraX command line bar before executing this script.")
            return

    from docx import Document
    summary_filename = "Batch_HBond_Summary.docx"
    doc = Document()
    
    # Document Title Setup
    doc.add_heading("Protein-Ligand Hydrogen Bond Summary Report", level=0)
    p_meta = doc.add_paragraph(f"Generated via ChimeraX Batch Automation Script\nCriteria: Geometric H-Bond Detection (distSlop: {global_slop} Å)")
    p_meta.runs[0].font.italic = True
    
    print(f"\n[INFO] Successfully initialized Word report layout: {summary_filename}")

    print("\n" + "="*50)
    print(" 2. EXECUTING BATCH PROCESS")
    print("="*50)

    # --- PHASE 2: LOOP THROUGH ALL FILES USING SAVED ANSWERS ---
    for filename in pdb_files:
        print("\n" + "*"*50)
        print(f" NOW PROCESSING: {filename}")
        print("*"*50)

        filename_base = os.path.splitext(filename)[0] 
        
        run(session, "close session")
        run(session, f"open {filename}")

        structures = session.models.list(type=AtomicStructure)
        if not structures:
            print(f"[ERROR] No structure loaded for {filename}. Skipping.")
            continue 

        # Apply saved chain choices
        current_file_chain_ids = []
        for s in structures:
            for c in s.chains:
                cid = c.chain_id.strip() if c.chain_id and c.chain_id.strip() else "NO-ID"
                if cid not in current_file_chain_ids:
                    current_file_chain_ids.append(cid)
                
        selected_chains = [c.strip() for c in global_chains_input.split(',')] if global_chains_input else current_file_chain_ids

        # Apply saved ligand choices
        current_file_ligands = []
        for s in structures:
            for r in s.residues:
                if r is not None and r.name not in std_amino_acids and r.name != "HOH":
                    current_file_ligands.append(r)
                    
        selected_ligands = []
        if global_ligands_input:
            try:
                indices = [int(x.strip()) for x in global_ligands_input.split(',')]
                for idx in indices:
                    if 0 <= idx < len(current_file_ligands):
                        selected_ligands.append(current_file_ligands[idx])
            except ValueError:
                pass

        # Build clean targets for proteins and ligands SEPARATELY
        chain_targets = []
        for c in selected_chains: 
            if c != "NO-ID": chain_targets.append(f"/{c}")
        protein_target = " | ".join(chain_targets)
        
        lig_targets = []
        for l in selected_ligands:
            c_id = l.chain_id.strip() if l.chain_id else "" 
            cid_part = f"/{c_id}" if c_id else ""
            lig_targets.append(f"#{l.structure.id_string}{cid_part}:{l.number}")
        ligand_target = " | ".join(lig_targets)

        unique_residues = set() 

        # --- RENDERING BLOCK ---
        try:
            run(session, "select clear")
            run(session, "hide atoms")
            run(session, "hide ribbons")
            run(session, "hide surfaces")
            
            if protein_target:
                run(session, f"show {protein_target} & protein ribbons")  
                run(session, f"color {protein_target} & protein cornflowerblue")
            
            if ligand_target:
                run(session, f"show {ligand_target} atoms")
                run(session, f"style {ligand_target} ball")
                run(session, f"color {ligand_target} yellow")
            
            run(session, "show ions")             
            run(session, "style ions sphere")
            run(session, "color ions green")
            
            if global_show_interact == 'y' and ligand_target:
                try:
                    # --- UPGRADE: Apply user-specified geometric cutoff tolerance ---
                    run(session, f"hbonds {ligand_target} restrict protein distSlop {global_slop}")
                    run(session, "color pbonds darkorange")
                    
                    # Dynamically extract strictly the residues involved in H-Bonds from the created pseudobonds
                    for m in session.models.list():
                        if hasattr(m, 'pseudobonds'):
                            for pb in m.pseudobonds:
                                a1, a2 = pb.atoms
                                if a1.residue.polymer_type == 1:
                                    unique_residues.add((a1.residue.chain_id, a1.residue.number))
                                if a2.residue.polymer_type == 1:
                                    unique_residues.add((a2.residue.chain_id, a2.residue.number))
                    
                    if unique_residues:
                        distinct_colors = ["crimson", "forestgreen", "darkmagenta", "darkorange", "teal", "firebrick", "indigo", "chocolate", "darkcyan", "purple", "olive", "indianred", "steelblue", "maroon", "darkblue"]
                        for idx, (chain_id, res_num) in enumerate(sorted(unique_residues)):
                            c_id = chain_id.strip() if chain_id else ""
                            res_spec = f"/{c_id}:{res_num}" if c_id else f":{res_num}"
                            assigned_color = distinct_colors[idx % len(distinct_colors)]
                            
                            run(session, f"color {res_spec} {assigned_color}")
                            run(session, f"label {res_spec} color {assigned_color} bgColor white height 1.2")
                            run(session, f"show {res_spec} atoms")
                            run(session, f"style {res_spec} stick")
                            
                    run(session, "select clear")

                    # --- APPEND TO WORD DOCX REPORT ---
                    doc.add_heading(f"Structure: {filename_base}", level=2)
                    
                    bonds_found = False
                    file_bonds = []
                    
                    for m in session.models.list():
                        if hasattr(m, 'pseudobonds'):
                            for pb in m.pseudobonds:
                                a1, a2 = pb.atoms
                                if a1.residue.polymer_type == 1:
                                    prot_atom, lig_atom = a1, a2
                                else:
                                    prot_atom, lig_atom = a2, a1
                                    
                                dist = math.dist(prot_atom.scene_coord, lig_atom.scene_coord)
                                prot_res_str = f"{prot_atom.residue.name} {prot_atom.residue.number}"
                                file_bonds.append((lig_atom.name, prot_res_str, prot_atom.name, f"{dist:.2f} Å"))
                                bonds_found = True
                    
                    if bonds_found:
                        # Construct native 4-column structured grid table
                        table = doc.add_table(rows=1, cols=4)
                        table.style = 'Table Grid'
                        
                        hdr_cells = table.rows[0].cells
                        hdr_cells[0].text = 'Ligand Atom'
                        hdr_cells[1].text = 'Protein Residue'
                        hdr_cells[2].text = 'Protein Atom'
                        hdr_cells[3].text = 'Distance (Å)'
                        
                        for lig_at, prot_res, prot_at, d_val in file_bonds:
                            row_cells = table.add_row().cells
                            row_cells[0].text = lig_at
                            row_cells[1].text = prot_res
                            row_cells[2].text = prot_at
                            row_cells[3].text = d_val
                    else:
                        p_none = doc.add_paragraph("No H-bonds detected within geometric criteria.")
                        p_none.runs[0].font.italic = True
                    
                    doc.add_paragraph("") # Spacer element

                except Exception as e:
                    print(f"-> No H-bonds mapped. Error: {e}")
                    doc.add_heading(f"Structure: {filename_base}", level=2)
                    p_err = doc.add_paragraph(f"Error evaluating interactions: {e}")
                    p_err.runs[0].font.italic = True
                
            run(session, "view")
            run(session, f"zoom {global_zoom}") 
            run(session, "lighting soft")
            run(session, "graphics silhouettes true")
            run(session, "set bgColor white")
            
            # 1. Save standard full-view image (Redirected inside Full_Renders folder)
            output_name_full = os.path.join(full_view_dir, f"{filename_base}_render_full.png")
            run(session, f'save "{output_name_full}" width 2400 height 2400 supersample 3')
            print(f"[SUCCESS] Saved full view to '{output_name_full}'")
            
            # 2. Save strict H-Bond zoomed view (Redirected inside Zoomed_Renders folder)
            if global_show_interact == 'y' and ligand_target:
                output_name_zoomed = os.path.join(zoomed_view_dir, f"{filename_base}_render_zoomed.png")
                
                # Dynamically frame the camera on ONLY the ligand and the strict H-bonded residues we found
                if unique_residues:
                    view_res_specs = [f"/{c_id.strip()}:{r_num}" if c_id.strip() else f":{r_num}" for c_id, r_num in unique_residues]
                    view_target = f"{ligand_target} | " + " | ".join(view_res_specs)
                    run(session, f"view {view_target}")
                else:
                    run(session, f"view {ligand_target}")
                
                run(session, "zoom 0.7") 
                run(session, f'save "{output_name_zoomed}" width 2400 height 2400 supersample 3')
                print(f"[SUCCESS] Saved zoomed view to '{output_name_zoomed}'")
            
        except Exception as e:
            print(f"[ERROR] Failed rendering {filename}: {e}")

    print("\n" + "="*50)
    print(" ALL FILES PROCESSED SUCCESSFULLY!")
    doc.save(summary_filename)
    print(f" Please check '{summary_filename}' for your Word Document data export.")
    print("="*50)
    run(session, "close session")

# Execute the function
process_batch(session)

```
</details>
