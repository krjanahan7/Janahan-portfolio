# 🐍 08_Python & R Scripts for Bioinformatics

![Tools](https://img.shields.io/badge/Tools-Biopython%20%7C%20pandas%20%7C%20ggplot2%20%7C%20tidyverse-blue)
![Language](https://img.shields.io/badge/Language-Python%20%7C%20R-orange)
![Skill](https://img.shields.io/badge/Skill-Data%20Analysis%20%7C%20Visualisation-yellow)

---

## 📌 Objective

This section collects reusable Python and R scripts I wrote for parsing biological file formats, wrangling result tables, and producing publication-ready plots.

---

## 🧠 Key Topics Covered

- FASTA/FASTQ/GenBank parsing with Biopython
- GC content, sequence statistics, and translation
- Data wrangling with pandas and dplyr
- Plotting with matplotlib, seaborn, and ggplot2
- Automating repetitive analysis steps
- Writing readable, reusable, documented functions

---

## 💻 Example Code

Python — sequence statistics:

    from Bio import SeqIO
    from Bio.SeqUtils import gc_fraction

    for record in SeqIO.parse("sequences.fasta", "fasta"):
        print(record.id, len(record.seq), round(gc_fraction(record.seq) * 100, 2))

Python — volcano plot:

    import pandas as pd, numpy as np, matplotlib.pyplot as plt
    df = pd.read_csv("DEG_results.csv")
    df["negLogP"] = -np.log10(df["padj"].clip(lower=1e-300))
    plt.scatter(df["log2FoldChange"], df["negLogP"], s=6, alpha=0.6)
    plt.xlabel("log2 Fold Change"); plt.ylabel("-log10 adjusted p")
    plt.savefig("volcano.png", dpi=300, bbox_inches="tight")

R — tidyverse summary and plot:

    library(tidyverse)
    deg <- read_csv("DEG_results.csv") %>%
      filter(padj < 0.05, abs(log2FoldChange) > 1)
    ggplot(deg, aes(log2FoldChange, -log10(padj))) +
      geom_point(alpha = 0.6) +
      theme_minimal() +
      labs(title = "Significant DEGs", x = "log2FC", y = "-log10 padj")

---

## ⚙️ Notes on Usage

| Item | Detail |
|---|---|
| Dependencies | `biopython`, `pandas`, `matplotlib`, `seaborn`; R: `tidyverse`, `DESeq2` |
| Input | Standard bioinformatics formats (FASTA, CSV, TSV, VCF) |
| Output | Cleaned tables and high-resolution figures (300 dpi) |
| Style | Functions documented with docstrings/comments for reuse |

---

## ⭐ Key Takeaway

Scripting turns a one-off analysis into a reproducible pipeline — clarity and reusability matter as much as correctness.

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

