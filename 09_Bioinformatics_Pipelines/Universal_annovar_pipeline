### 📌 Universal NGS Variant Calling & Annotation Pipeline (ANNOVAR)

[![View Script](https://img.shields.io/badge/VIEW_SCRIPT-UNIVERSAL__ANNOVAR__PIPELINE.PY-blue?style=for-the-badge&logo=python)](./universal_annovar_pipeline.py)
[![Status](https://img.shields.io/badge/STATUS-EXECUTABLE-brightgreen?style=for-the-badge)](#)

<details>
<summary>🔍 <b>Click to View Script Details & Usage Header</b></summary>

```python
# ==============================================================================
# SCRIPT NAME : universal_annovar_pipeline.py
# AUTHOR      : K. R. Janahan
# DESCRIPTION : Automated NGS pipeline for paired-end Illumina read alignment,
#               BAM sorting/indexing, FreeBayes variant calling, and ANNOVAR annotation.
# DEPENDENCIES: BWA, Samtools, FreeBayes, Perl, ANNOVAR
# USAGE       :
# conda activate universal_annovar_pipeline
# python3 universal_annovar_pipeline.py \
#   --r1 /path/to/new_sample_R1.fastq.gz \
#   --r2 /path/to/new_sample_R2.fastq.gz \
#   --sample NewSampleName \
#   --outdir ./results_new_sample \
#   --reference /path/to/your_genome.fna \
#   --annovar_dir /home/janahan/ANNOVAR/ \
#   --annovar_db /home/janahan/ANNOVAR/humandb/ \
#   --buildver hg38
# ==============================================================================
```
</details>

<details>
<summary>💻 <b>Click to View / Copy Full Source Code</b></summary>

```python

import argparse
import subprocess
import os
import sys
import shutil

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

def run_command(command, step_name):
    """Executes a shell command and handles errors robustly."""
    print(f"\n[INFO] Starting: {step_name}")
    print(f"[CMD] {' '.join(command)}")
    try:
        subprocess.run(command, check=True, text=True, capture_output=True)
        print(f"[SUCCESS] {step_name} completed.")
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] {step_name} failed with exit code {e.returncode}.")
        print(f"[STDERR]\n{e.stderr}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Universal NGS Pipeline")
    
    # Input/Output Arguments
    parser.add_argument("--r1", required=True, help="Path to Read 1 FASTQ (compressed or uncompressed)")
    parser.add_argument("--r2", required=True, help="Path to Read 2 FASTQ (compressed or uncompressed)")
    parser.add_argument("--sample", required=True, help="Sample Name")
    parser.add_argument("--outdir", required=True, help="Output directory path")
    
    # Resource Arguments
    parser.add_argument("--reference", required=True, help="Path to reference genome (.fna or .fna.gz)")
    parser.add_argument("--annovar_dir", required=True, help="Path to Annovar installation")
    parser.add_argument("--annovar_db", required=True, help="Path to Annovar humandb")
    parser.add_argument("--buildver", default="hg38", help="Genome build version")
    parser.add_argument("--threads", default="16", help="Number of threads to use")
    
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # --- STEP 1: REFERENCE GENOME HANDLING (Decompression & Renaming) ---
    print("\n[INFO] --- STEP 1: Reference Genome Pre-processing ---")
    standardized_ref = os.path.join(args.outdir, "ref.fna")

    if args.reference.endswith(".gz"):
        print(f"[INFO] Reference is compressed (.gz). Decompressing and renaming directly to: {standardized_ref}")
        with open(standardized_ref, "w") as out_ref:
            subprocess.run(["gunzip", "-c", args.reference], stdout=out_ref, check=True)
    else:
        print(f"[INFO] Reference is already uncompressed. Standardizing name to: {standardized_ref}")
        shutil.copy(args.reference, standardized_ref)

    # Mandatory Reference Indexing
    run_command(["bwa", "index", standardized_ref], "BWA Reference Indexing")
    run_command(["samtools", "faidx", standardized_ref], "Samtools Reference Indexing")

    # --- STEP 2: FASTQ DECOMPRESSION ---
    print("\n[INFO] --- STEP 2: FastQ Reads Pre-processing ---")
    unzipped_r1 = os.path.join(args.outdir, f"{args.sample}_R1_uncompressed.fastq")
    unzipped_r2 = os.path.join(args.outdir, f"{args.sample}_R2_uncompressed.fastq")
    
    print(f"[INFO] Decompressing R1 read files...")
    with open(unzipped_r1, "w") as out_r1:
        subprocess.run(["gunzip", "-c", args.r1], stdout=out_r1, check=True)
        
    print(f"[INFO] Decompressing R2 read files...")
    with open(unzipped_r2, "w") as out_r2:
        subprocess.run(["gunzip", "-c", args.r2], stdout=out_r2, check=True)

    # Define downstream intermediate and final file paths
    sam_file = os.path.join(args.outdir, f"{args.sample}.sam")
    bam_file = os.path.join(args.outdir, f"{args.sample}.bam")
    sorted_bam = os.path.join(args.outdir, f"{args.sample}.sorted.bam")
    vcf_file = os.path.join(args.outdir, f"{args.sample}.vcf")
    annovar_out_prefix = os.path.join(args.outdir, f"{args.sample}_annotated")

    # --- STEP 3: ALIGNMENT (BWA MEM) ---
    bwa_cmd = [
        "bwa", "mem", "-t", args.threads, 
        "-R", f"@RG\\tID:{args.sample}\\tSM:{args.sample}\\tPL:ILLUMINA",
        standardized_ref, unzipped_r1, unzipped_r2
    ]
    
    print(f"\n[INFO] --- STEP 3: BWA Alignment ---")
    with open(sam_file, "w") as out_sam:
        subprocess.run(bwa_cmd, stdout=out_sam, check=True)
    print(f"[SUCCESS] BWA Alignment completed.")

    # --- STEP 4: SAM TO BAM & SORTING (Samtools) ---
    print(f"\n[INFO] --- STEP 4: Samtools Processing ---")
    samtools_view = ["samtools", "view", "-@", args.threads, "-bS", sam_file]
    with open(bam_file, "w") as out_bam:
         subprocess.run(samtools_view, stdout=out_bam, check=True)

    samtools_sort = ["samtools", "sort", "-@", args.threads, "-o", sorted_bam, bam_file]
    run_command(samtools_sort, "Samtools Sort")

    samtools_index = ["samtools", "index", sorted_bam]
    run_command(samtools_index, "Samtools Index")

    # Clean up massive intermediate text/fastq files to save disk storage
    os.remove(sam_file)
    os.remove(bam_file) # Removing unsorted BAM to save space
    os.remove(unzipped_r1)
    os.remove(unzipped_r2)

    # --- STEP 5: VARIANT CALLING (FreeBayes) ---
    print(f"\n[INFO] --- STEP 5: Freebayes Variant Calling ---")
    freebayes_cmd = ["freebayes", "-f", standardized_ref, sorted_bam]
    with open(vcf_file, "w") as out_vcf:
        subprocess.run(freebayes_cmd, stdout=out_vcf, check=True)
    print(f"[SUCCESS] Freebayes Variant Calling completed.")

    # --- STEP 6: ANNOVAR ANNOTATION ---
    print(f"\n[INFO] --- STEP 6: Annovar Functional Annotation ---")
    table_annovar_script = os.path.join(args.annovar_dir, "table_annovar.pl")
    protocols = "refGene,clinvar_20221231,gnomad211_exome" 
    operations = "g,f,f"

    annovar_cmd = [
        "perl", table_annovar_script, vcf_file, args.annovar_db,
        "-buildver", args.buildver,
        "-out", annovar_out_prefix,
        "-remove", 
        "-protocol", protocols,
        "-operation", operations,
        "-nastring", ".", 
        "-vcfinput" 
    ]
    
    run_command(annovar_cmd, "Annovar Annotation")

    print(f"\n[PIPELINE COMPLETE] Your final annotated variant files are ready at: {annovar_out_prefix}.*")

if __name__ == "__main__":
    main()

```
</details>
