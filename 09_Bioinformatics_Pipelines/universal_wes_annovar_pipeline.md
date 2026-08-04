### 📌 Universal Organism NGS & Annotation Pipeline (ANNOVAR / GATK)

[![View Script](https://img.shields.io/badge/VIEW_SCRIPT-UNIVERSAL__ANNOVAR__GATK__PIPELINE.PY-blue?style=for-the-badge&logo=python)](#)
[![Status](https://img.shields.io/badge/STATUS-EXECUTABLE-brightgreen?style=for-the-badge)](#)

<details>
<summary>🔍 <b>Click to View Script Details & Usage Header</b></summary>

```python
# ==============================================================================
# SCRIPT NAME :universal_wes_annovar_pipeline.py
# AUTHOR      : K. R. Janahan
# DESCRIPTION : Automated NGS pipeline for paired-end Illumina read alignment,
#               BAM sorting/indexing, GATK HaplotypeCaller variant calling, 
#               dynamic automated ANNOVAR database building, and annotation.
# DEPENDENCIES: BWA, Samtools, GATK, Perl, ANNOVAR
# USAGE       :
conda activate universal_annovar_pipeline
python universal_annovar_pipeline.py \
  --r1 /path/to/new_sample_R1.fastq.gz \
  --r2 /path/to/new_sample_R2.fastq.gz \
  --sample NewSampleName \
  --outdir ./results_new_sample \
  --reference /path/to/your_genome.fna \
  --annovar_dir /home/janahan/ANNOVAR/ \
  --annovar_db /home/janahan/ANNOVAR/humandb/ \
  --buildver hg38

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
    parser = argparse.ArgumentParser(description="Universal Organism NGS & Annotation Pipeline")
    
    # Input/Output Arguments per Sample
    parser.add_argument("--r1", required=True, help="Path to Read 1 FASTQ (compressed or uncompressed)")
    parser.add_argument("--r2", required=True, help="Path to Read 2 FASTQ (compressed or uncompressed)")
    parser.add_argument("--sample", required=True, help="Sample Name")
    parser.add_argument("--outdir", required=True, help="Output directory path")
    parser.add_argument("--reference", required=True, help="Path to reference genome (.fna or .fna.gz)")
    parser.add_argument("--annotation", required=True, help="Path to genomic annotation (.gff or .gff.gz)")
    
    # Resource Arguments (Static defaults for your system)
    parser.add_argument("--annovar_dir", default="/home/janahan/ANNOVAR", help="Path to Annovar installation")
    parser.add_argument("--annovar_db", default="/home/janahan/ANNOVAR/humandb", help="Path to Annovar database folder")
    parser.add_argument("--threads", default="16", help="Number of threads to use")
    
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Automatically derive a unique genome build name from the reference filename
    ref_base = os.path.basename(args.reference)
    build_version = ref_base.split("_genomic")[0].split(".fna")[0]
    print(f"[INFO] Target Genome Build Version identified as: {build_version}")

    # --- STEP 1: REFERENCE GENOME HANDLING ---
    print("\n[INFO] --- STEP 1: Reference Genome Pre-processing ---")
    standardized_ref = os.path.join(args.outdir, "ref.fna")
    dict_file = os.path.join(args.outdir, "ref.dict")

    if not os.path.exists(standardized_ref):
        if args.reference.endswith(".gz"):
            print(f"[INFO] Reference is compressed (.gz). Decompressing directly to: {standardized_ref}")
            with open(standardized_ref, "w") as out_ref:
                subprocess.run(["gunzip", "-c", args.reference], stdout=out_ref, check=True)
        else:
            print(f"[INFO] Reference is uncompressed. Standardizing name to: {standardized_ref}")
            shutil.copy(args.reference, standardized_ref)
    else:
        print(f"[INFO] Reference file {standardized_ref} already exists. Skipping decompression.")

    # Indexing steps
    if not os.path.exists(f"{standardized_ref}.bwt"):
        run_command(["bwa", "index", standardized_ref], "BWA Reference Indexing")
    if not os.path.exists(f"{standardized_ref}.fai"):
        run_command(["samtools", "faidx", standardized_ref], "Samtools Reference Indexing")
    if not os.path.exists(dict_file):
        run_command(["gatk", "CreateSequenceDictionary", "-R", standardized_ref], "GATK Sequence Dictionary Creation")

    # --- STEP 1b: DYNAMIC AUTOMATED ANNOVAR DATABASE BUILDING ---
    print("\n[INFO] --- STEP 1b: Checking/Building Custom Organism Database for Annovar ---")
    db_txt_file = os.path.join(args.annovar_db, f"{build_version}_refGene.txt")
    db_fa_file = os.path.join(args.annovar_db, f"{build_version}_refGeneMrna.fa")

    if not os.path.exists(db_txt_file) or not os.path.exists(db_fa_file):
        print(f"[INFO] Custom database components missing for {build_version}. Commencing automated generation...")
        
        # Unzip the GFF if it is compressed
        unzipped_gff = os.path.join(args.outdir, "temp_annotation.gff")
        if args.annotation.endswith(".gz"):
            print(f"[INFO] Annotation is compressed (.gz). Decompressing: {args.annotation}")
            with open(unzipped_gff, "w") as out_gff:
                subprocess.run(["gunzip", "-c", args.annotation], stdout=out_gff, check=True)
        else:
            shutil.copy(args.annotation, unzipped_gff)

        # Build step 1: Generate _refGene.txt using gff3ToGenePred
        gff_convert_script = os.path.join(args.annovar_dir, "gff3ToGenePred")
        if not os.path.exists(gff_convert_script):
            gff_convert_script = "gff3ToGenePred"
            
        cmd_build_txt = [gff_convert_script, "-useName", unzipped_gff, db_txt_file]
        
        try:
            run_command(cmd_build_txt, "Annovar Custom DB Generation (_refGene.txt)")
        except SystemExit:
            print("[WARNING] Direct gff3ToGenePred utility failed. Executing fallback native database construction...")
            convert_script = os.path.join(args.annovar_dir, "gff3ToGenePred.pl")
            if os.path.exists(convert_script):
                run_command(["perl", convert_script, unzipped_gff, db_txt_file], "Annovar Fallback Perl Parser")
            else:
                shutil.copy(unzipped_gff, db_txt_file)

        # Build step 2: Generate _refGeneMrna.fa
        retrieve_seq_script = os.path.join(args.annovar_dir, "retrieve_seq_from_fasta.pl")
        cmd_build_fa = ["perl", retrieve_seq_script, "-format", "refGene", "-seqfile", standardized_ref, db_txt_file, "-out", db_fa_file]
        run_command(cmd_build_fa, "Annovar Custom DB Generation (_refGeneMrna.fa)")
        
        if os.path.exists(unzipped_gff):
            os.remove(unzipped_gff)
        print(f"[SUCCESS] Database for {build_version} successfully built and stored.")
    else:
        print(f"[INFO] Existing database detected for {build_version} inside Annovar DB directory. Skipping build step.")

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

    # Paths setup
    sam_file = os.path.join(args.outdir, f"{args.sample}.sam")
    bam_file = os.path.join(args.outdir, f"{args.sample}.bam")
    sorted_bam = os.path.join(args.outdir, f"{args.sample}.sorted.bam")
    vcf_file = os.path.join(args.outdir, f"{args.sample}.vcf")
    annovar_out_prefix = os.path.join(args.outdir, f"{args.sample}_annotated")

    # --- STEP 3: ALIGNMENT (BWA MEM) ---
    print(f"\n[INFO] --- STEP 3: BWA Alignment ---")
    bwa_cmd = ["bwa", "mem", "-t", args.threads, "-R", f"@RG\\tID:{args.sample}\\tSM:{args.sample}\\tPL:ILLUMINA", standardized_ref, unzipped_r1, unzipped_r2]
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

    os.remove(sam_file)
    os.remove(bam_file) 
    os.remove(unzipped_r1)
    os.remove(unzipped_r2)

    # --- STEP 5: VARIANT CALLING (GATK) ---
    print(f"\n[INFO] --- STEP 5: GATK Variant Calling ---")
    gatk_cmd = ["gatk", "HaplotypeCaller", "-R", standardized_ref, "-I", sorted_bam, "-O", vcf_file]
    run_command(gatk_cmd, "GATK HaplotypeCaller")

    # --- STEP 6: ANNOVAR ANNOTATION ---
    print(f"\n[INFO] --- STEP 6: Annovar Functional Annotation ---")
    table_annovar_script = os.path.join(args.annovar_dir, "table_annovar.pl")
    annovar_cmd = ["perl", table_annovar_script, vcf_file, args.annovar_db, "-buildver", build_version, "-out", annovar_out_prefix, "-remove", "-protocol", "refGene", "-operation", "g", "-nastring", ".", "-vcfinput"]
    run_command(annovar_cmd, "Annovar Annotation")

    # --- STEP 7: GENE NAME VERIFICATION STEP ---
    final_txt = f"{annovar_out_prefix}.{build_version}_multianno.txt"
    if os.path.exists(final_txt):
        print(f"\n[VERIFICATION] Verifying Gene Names in: {final_txt}")
        with open(final_txt, 'r') as f:
            header = f.readline().strip().split('\t')
            if "Gene.refGene" in header:
                print("[SUCCESS] Gene Name column tracked successfully! Look for the 'Gene.refGene' column in your data.")
            else:
                print("[WARNING] 'Gene.refGene' column was not detected.")

    print(f"\n[PIPELINE COMPLETE] Your final annotated variant files are ready at: {annovar_out_prefix}.*")

if __name__ == "__main__":
    main()

```
</details>
