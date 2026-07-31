# 💊 07_Molecular Docking & Drug Discovery

![Tools](https://img.shields.io/badge/Tools-AutoDock%20Vina%20%7C%20Open%20Babel%20%7C%20PyMOL-blue)
![Language](https://img.shields.io/badge/Language-Bash%20%7C%20Python-orange)
![Stage](https://img.shields.io/badge/Stage-Computer--Aided%20Drug%20Design-yellow)

---

## 📌 Objective

This section documents structure-based virtual screening: preparing receptor and ligand files, defining a binding site, running docking simulations, and interpreting binding affinities and interactions.

---

## 🧠 Key Topics Covered

- Receptor preparation (water removal, hydrogens, charges)
- Ligand preparation and energy minimisation
- Grid box definition around the active site
- Docking with AutoDock Vina
- Binding affinity and pose ranking
- Protein-ligand interaction analysis

---

## 💻 Commands Used

Prepare files:

    obabel ligand.sdf -O ligand.pdbqt --gen3d -p 7.4
    obabel receptor.pdb -O receptor.pdbqt -xr

Configuration file (`config.txt`):

    receptor = receptor.pdbqt
    ligand = ligand.pdbqt
    center_x = 12.5
    center_y = -4.8
    center_z = 30.2
    size_x = 24
    size_y = 24
    size_z = 24
    exhaustiveness = 16
    num_modes = 9

Run docking:

    vina --config config.txt --out docked_poses.pdbqt --log docking_log.txt
    obabel docked_poses.pdbqt -O docked_poses.pdb -m

---

## ⚙️ Parameter Explanations

| Parameter | Meaning |
|---|---|
| `--gen3d` | Generates 3D coordinates for the ligand |
| `-p 7.4` | Adds hydrogens at physiological pH |
| `-xr` | Treats the receptor as rigid |
| `center_x/y/z` | Coordinates of the binding-site centre |
| `size_x/y/z` | Grid box dimensions in Ångströms |
| `exhaustiveness` | Search thoroughness — higher is slower but more reliable |
| `num_modes` | Number of binding poses to report |

---

## 📊 Results & Interpretation

- Binding affinity is reported in kcal/mol; more negative values indicate stronger predicted binding.
- Poses within 2 kcal/mol of the best score are inspected for consistent interactions.
- Hydrogen bonds, hydrophobic contacts, and pi-stacking with key residues are recorded.

---

## ⭐ Key Takeaway

Docking scores rank hypotheses, not truths — the credibility of a pose comes from chemically sensible interactions with known active-site residues.

---

## 📬 Contact & Links

* **GitHub:** [krjanahan7](https://github.com/krjanahan7)
* **LinkedIn:** [Janahan KR](https://www.linkedin.com/in/janahan-kr-198790375)
* **Email:** [krjanahan7@gmail.com](mailto:krjanahan7@gmail.com)

---

## 🔗 Back to Main Portfolio

📂 **Portfolio Home:** [Janahan-Portfolio](https://github.com/Janahan10/Janahan-Portfolio)

---

*Generated for professional bioinformatics portfolio presentation.*

