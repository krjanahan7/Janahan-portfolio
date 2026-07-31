# 🧬 FastQC & Quality Control

![Tool](https://img.shields.io/badge/Tools-FastQC%20%7C%20MultiQC-blue)
![Language](https://img.shields.io/badge/Language-Bash-orange)
![Pipeline](https://img.shields.io/badge/Stage-Pre--processing-yellow)

---

## 📌 Objective

This section documents my workflow for assessing raw sequencing data quality using **FastQC** and **MultiQC**. Ensuring high read quality is a mandatory first step before any downstream alignment, assembly, variant calling, or differential expression analysis. The goal is to identify technical artifacts such as adapter contamination, low-quality bases, and overrepresented sequences, and to document any quality control decisions made before proceeding with the pipeline.

---

## 🛠️ Key Topics Covered

- **Per-Base Sequence Quality:** Interpreting Phred quality scores across all read positions.
- **Per-Sequence Quality Scores:** Detecting reads with universally low quality.
- **GC Content & Sequence Duplication Levels:** Checking for contamination, library preparation bias, or PCR duplicates.
- **Adapter Contamination:** Identifying overrepresented sequences and adapter artifacts.
- **Read Length Distribution:** Validating expected fragment sizes and trimming needs.
- **Multi-Report Aggregation:** Using MultiQC to summarize quality reports across dozens of samples simultaneously.
- **QC-Based Trimming Decisions:** Deciding when to trim, filter, or remove low-quality reads.

---

## 💻 Workflow / Commands

### Step 1: Create Output Directory for Reports

    mkdir -p 02_FastQC/results/fastqc_reports
    mkdir -p 02_FastQC/results/multiqc_summary

### Step 2: Run FastQC on a Single FASTQ File

    fastqc -o ./results/fastqc_reports/ -t 4 sample_R1.fastq.gz

### Step 3: Run FastQC in Batch on All FASTQ Files

    for file in ./data/*.fastq.gz
    do
        echo "Running FastQC on $file..."
        fastqc -o ./results/fastqc_reports/ -t 4 "$file"
    done

### Step 4: Aggregate All Reports with MultiQC

    multiqc ./results/fastqc_reports/ -o ./results/multiqc_summary/

### Step 5: Optional — Trim Low-Quality Bases and Adapters with Trimmomatic

    trimmomatic PE -threads 4 \
        sample_R1.fastq.gz sample_R2.fastq.gz \
        sample_R1_paired.fq.gz sample_R1_unpaired.fq.gz \
        sample_R2_paired.fq.gz sample_R2_unpaired.fq.gz \
        ILLUMINACLIP: TruSeq3-PE.fa:2:30:10 \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36

---

## 🔑 Key Parameters & Explanation

| Parameter / Command | Description |
|---------------------|-------------|
| `fastqc` | Tool for high-throughput sequence QC analysis and HTML report generation. |
| `-o` | Output directory where FastQC report files (`.html` and `.zip`) are saved. |
| `-t 4` | Number of threads allocated to FastQC for parallel processing. |
| `*.fastq.gz` | Wildcard pattern to process all compressed FASTQ files in a directory. |
| `multiqc` | Tool that searches a directory for QC reports and generates a single interactive summary. |
| `Trimmomatic PE` | Paired-end trimming tool to remove adapters and low-quality bases. |
| `ILLUMINACLIP` | Removes adapter sequences using a provided adapter FASTA file. |
| `LEADING:20` | Removes low-quality bases from the start of reads below Q20. |
| `TRAILING:20` | Removes low-quality bases from the end of reads below Q20. |
| `SLIDINGWINDOW:4:20` | Scans a 4-base window and cuts when average quality drops below Q20. |
| `MINLEN:36` | Drops reads shorter than 36 bases after trimming. |

> **Quality Threshold:** A Phred quality score (`Q`) above 30 means a base call accuracy of ≥ 99.9%, which is generally considered acceptable for downstream analysis.

---

## 📊 Sample Output / Results

This section should include screenshots and summary files from the actual QC run. Examples to include:

- **FastQC Per-Base Sequence Quality Plot:** Screenshot showing quality scores across all read positions.
- **FastQC Per-Sequence GC Content Plot:** Distribution of GC content per read compared to theoretical distribution.
- **FastQC Adapter Content Plot:** Evidence of adapter contamination across read lengths.
- **MultiQC Summary HTML:** Aggregated report across all samples with pass/warn/fail flags.

### Example Results Summary Table

| Sample ID | Total Reads | Mean Q Score | GC % | Adapter Content | Status |
|-----------|-------------|--------------|------|-----------------|--------|
| sample_R1 | 25,000,000 | 34.5 | 49.2% | None | ✅ Pass |
| sample_R2 | 25,000,000 | 33.8 | 49.0% | Low | ⚠️ Trim |

---

## 📁 Files in This Folder

| File / Directory | Description |
|------------------|-------------|
| `scripts/` | Reusable shell scripts for batch FastQC and MultiQC execution. |
| `data/` | Example raw FASTQ files or sample lists used for QC. |
| `results/fastqc_reports/` | Individual FastQC HTML and `.zip` reports. |
| `results/multiqc_summary/` | Aggregated MultiQC HTML report and data. |
| `figures/` | Screenshots of key plots used in the portfolio. |
| `README.md` | This documentation file. |

---

## ⭐ Key Takeaway

High-quality raw reads are the foundation of every reliable bioinformatics analysis. A rigorous FastQC/MultiQC pre-processing step prevents garbage-in-garbage-out scenarios and ensures downstream alignment, quantification, and variant calling results are biologically meaningful.

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
