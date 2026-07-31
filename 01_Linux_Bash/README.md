# 🐧 01_Linux & Bash for Bioinformatics

![Tool](https://img.shields.io/badge/Tools-Linux%20%7C%20Bash%20%7C%20Shell%20Scripting-blue)
![Language](https://img.shields.io/badge/Language-Bash%2FShell-orange)
![Skill](https://img.shields.io/badge/Skill-Command%20Line%20%7C%20Automation-yellow)

---

## 📌 Objective

This section documents my practical experience using the **Linux command line** and **Bash scripting** for file management, data processing, and automation in bioinformatics workflows. The goal is to demonstrate proficiency in navigating Unix-based systems, manipulating large biological datasets, and writing reusable shell scripts that streamline repetitive analysis tasks.

Mastering Linux and Bash is essential for working in computational biology, as most high-performance computing (HPC) clusters, genomics pipelines, and open-source bioinformatics tools run in Linux environments.

---

## 🛠️ Key Topics Covered

- **File Navigation & Management:** Using `ls`, `cd`, `mkdir`, `cp`, `mv`, `rm`, and `find` to organize data.
- **Text Manipulation & Searching:** Using `grep`, `awk`, `sed`, `cut`, `sort`, `uniq`, and `wc` to extract and summarize data.
- **Sequence Data Inspection:** Viewing FASTQ/FASTA headers and counts using `head`, `tail`, `wc -l`, and `zcat`.
- **Permission & Environment Control:** Understanding file permissions (`chmod`, `chown`) and managing Conda/Mamba environments.
- **Pipeline Automation:** Writing `for` loops, `while` loops, and shell scripts to process multiple samples simultaneously.
- **Redirect & Piping:** Combining commands with `|`, `>`, `>>`, and `2>` for efficient data streams.
- **Job Scheduling Basics:** Preparing SLURM-ready scripts for running analyses on HPC clusters.

---

## 💻 Workflow / Commands

### Step 1: Navigate and Inspect Project Files

    # List all files with details, including hidden files
    ls -lha
    
    # Check current directory
    pwd
    
    # Create a new project folder for analysis
    mkdir -p 01_Linux_Bash/data 01_Linux_Bash/scripts 01_Linux_Bash/results

### Step 2: Count Total Reads in a FASTQ File

    # Each FASTQ record has 4 lines, so divide total lines by 4
    echo $(( $(wc -l < data/sample.fastq) / 4 ))

### Step 3: Extract Header Lines from a FASTA File

    grep "^>" data/genome.fasta > results/headers.txt

### Step 4: Count Sequences in a FASTA File

    grep -c "^>" data/genome.fasta

### Step 5: Compress and Decompress Large Files

    # Compress a FASTQ file
    gzip data/sample.fastq
    
    # View the first few lines of a compressed FASTQ without extracting
    zcat data/sample.fastq.gz | head -n 12

### Step 6: Batch Process Multiple FASTQ Files with a For Loop

    for file in data/*.fastq.gz
    do
        echo "Processing $file..."
        base=$(basename "$file" .fastq.gz)
        echo $base >> results/processed_samples.txt
    done

### Step 7: Extract Specific Columns with awk and cut

    # Extract the first and third columns from a tab-delimited file
    awk -F'\t' '{print $1, $3}' data/annotations.tsv > results/gene_ids.txt
    
    # Extract the sequence ID column from a FASTA header
    grep "^>" data/genome.fasta | cut -d' ' -f1 > results/seq_ids.txt

### Step 8: Search and Replace Patterns with sed

    # Replace all occurrences of a sample name in a file
    sed 's/old_sample/new_sample/g' data/sample_list.txt > results/sample_list_cleaned.txt

### Step 9: Make a Script Executable and Run It

    chmod +x scripts/run_analysis.sh
    ./scripts/run_analysis.sh

### Step 10: Sample SLURM Job Script for HPC Clusters

    #!/bin/bash
    #SBATCH --job-name=linux_bioinfo
    #SBATCH --output=results/job_%j.out
    #SBATCH --error=results/job_%j.err
    #SBATCH --time=01:00:00
    #SBATCH --cpus-per-task=4
    #SBATCH --mem=8G
    
    echo "Job started on $(date)"
    bash scripts/run_analysis.sh
    echo "Job finished on $(date)"

---

## 🔑 Key Parameters & Explanation

| Parameter / Command | Description |
|---------------------|-------------|
| `ls -lha` | Lists files in long format with human-readable sizes and hidden files. |
| `mkdir -p` | Creates directories recursively without error if they already exist. |
| `wc -l < file` | Counts the number of lines in a file (input redirection). |
| `grep "^>"` | Matches lines starting with `>`, used to extract FASTA headers. |
| `grep -c` | Counts the number of matching lines. |
| `zcat` | Displays the contents of a gzip-compressed file without extracting it. |
| `basename` | Strips directory and optionally suffix from a filename. |
| `awk -F'\t' '{print $1, $3}'` | Splits lines by tab and prints selected columns. |
| `cut -d' ' -f1` | Splits lines by a delimiter and prints the first field. |
| `sed 's/old/new/g'` | Substitutes all occurrences of a pattern in a file. |
| `chmod +x` | Makes a script file executable. |
| `|` | Pipe operator that passes the output of one command as input to another. |
| `>` / `>>` | Redirects output to a file (overwrite / append). |
| `2>` | Redirects standard error messages to a file. |

---

## 📊 Sample Output / Results

This section should include examples of command outputs and scripts in action. Examples to include:

- **Terminal Screenshot:** Showing a sequence of successful Linux commands.
- **Sample Count Summary:** Output table showing read counts per FASTQ file.
- **Script Output:** Log from a batch processing script.
- **FASTA Header Preview:** Screenshot of extracted headers from a genome file.

### Example Results Summary Table

| File | Total Reads | Total Bases | Command Used |
|------|-------------|-------------|--------------|
| sample_1.fastq.gz | 12,500,000 | 1,250,000,000 | `zcat file | wc -l` / 4 |
| sample_2.fastq.gz | 14,200,000 | 1,420,000,000 | `zcat file | wc -l` / 4 |
| genome.fasta | 8,450 sequences | 3,200,000,000 bp | `grep -c "^>"` |

---

## 📁 Files in This Folder

| File / Directory | Description |
|------------------|-------------|
| `scripts/` | Reusable shell scripts for automation and batch processing. |
| `data/` | Example FASTQ, FASTA, and TSV files used for practice. |
| `results/` | Output files generated during analysis. |
| `notes/` | Reference notes on Linux commands and bash syntax. |
| `README.md` | This documentation file. |

---

## ⭐ Key Takeaway

Linux and Bash proficiency is the backbone of bioinformatics. Efficient command-line skills enable reproducible data handling, automation of multi-sample workflows, and seamless integration with HPC clusters and workflow management systems.

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
