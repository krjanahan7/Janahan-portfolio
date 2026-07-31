# 🧱 06_Structural Biology & Protein Modelling

![Tools](https://img.shields.io/badge/Tools-PyMOL%20%7C%20SWISS--MODEL%20%7C%20AlphaFold%20%7C%20Ramachandran-blue)
![Language](https://img.shields.io/badge/Language-Python%20%7C%20PyMOL%20Script-orange)
![Stage](https://img.shields.io/badge/Stage-Structural%20Analysis-yellow)

---

## 📌 Objective

This section documents protein structure retrieval, homology modelling, structure validation, and visualisation for structure-function interpretation.

---

## 🧠 Key Topics Covered

- PDB structure retrieval and file format understanding
- Homology modelling and template selection
- AlphaFold predicted structures and pLDDT confidence
- Model validation (Ramachandran plot, QMEAN)
- Structure superposition and RMSD
- Publication-quality visualisation

---

## 💻 Commands Used

Fetch and inspect:

    wget https://files.rcsb.org/download/1ABC.pdb

PyMOL session:

    fetch 1ABC
    hide everything
    show cartoon
    color skyblue, chain A
    select active_site, resi 45+78+120
    show sticks, active_site
    align model.pdb, 1ABC
    png structure_figure.png, width=1600, height=1200, dpi=300, ray=1

Biopython parsing:

    from Bio.PDB import PDBParser
    structure = PDBParser(QUIET=True).get_structure("prot", "1ABC.pdb")
    for chain in structure[0]:
        print(chain.id, len(list(chain.get_residues())))

---

## ⚙️ Parameter Explanations

| Parameter | Meaning |
|---|---|
| `fetch` | Downloads structure directly from the PDB |
| `align` | Superposes two structures and reports RMSD |
| `ray=1` | Ray-traced rendering for high-quality figures |
| `pLDDT` | AlphaFold per-residue confidence score (0–100) |
| `RMSD` | Root-mean-square deviation between aligned atoms |

---

## 📊 Results & Interpretation

- RMSD below ~2 Å indicates strong structural similarity to the template.
- Ramachandran analysis confirms most residues lie in favoured regions.
- Active-site residues are mapped and visualised for downstream docking.

---

## ⭐ Key Takeaway

A protein model is a hypothesis — validation metrics decide whether it can support biological conclusions.

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
