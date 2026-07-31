# 🔬 05_RNA-Seq & Differential Expression

![Tools](https://img.shields.io/badge/Tools-HISAT2%20%7C%20featureCounts%20%7C%20DESeq2-blue)
![Language](https://img.shields.io/badge/Language-R%20%7C%20Bash-orange)
![Stage](https://img.shields.io/badge/Stage-Transcriptomics-yellow)

---

## 📌 Objective

This section documents a complete RNA-Seq workflow, from spliced alignment and gene-level quantification to statistical detection of differentially expressed genes and biological interpretation.

---

## 🧠 Key Topics Covered

- Spliced alignment and transcript quantification
- Count matrix generation
- Normalisation and dispersion estimation
- Differential expression testing with DESeq2
- PCA, MA plots, volcano plots, heatmaps
- Functional enrichment (GO / KEGG)

---

## 💻 Commands Used

Quantification:

    featureCounts -T 4 -p -a annotation.gtf -o counts.txt *.sorted.bam

Pseudo-alignment alternative:

    salmon quant -i salmon_index -l A -1 R1.fq.gz -2 R2.fq.gz -p 4 -o quants/sample1

DESeq2 in R:

    library(DESeq2)
    counts <- read.delim("counts.txt", row.names = 1, comment.char = "#")
    coldata <- read.csv("metadata.csv", row.names = 1)
    dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", "treated", "control"))
    write.csv(as.data.frame(res[order(res$padj), ]), "DEG_results.csv")

---

## ⚙️ Parameter Explanations

| Parameter | Meaning |
|---|---|
| `-p` | Count fragments instead of reads for paired-end data |
| `-a annotation.gtf` | Gene model used for feature assignment |
| `design = ~ condition` | Statistical model comparing experimental groups |
| `padj` | Benjamini-Hochberg adjusted p-value (FDR) |
| `log2FoldChange` | Effect size of expression change |
| `rowSums >= 10` | Pre-filter removes low-count noise genes |

---

## 📊 Results & Interpretation

- PCA separates treated and control samples along PC1, confirming a biological effect.
- Significant genes are defined as `padj < 0.05` and `|log2FC| > 1`.
- Enrichment analysis links the DEG list back to pathways and biological processes.

---

## ⭐ Key Takeaway

Differential expression is a statistical claim, not a raw count comparison — replication and multiple-testing correction are essential.

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

