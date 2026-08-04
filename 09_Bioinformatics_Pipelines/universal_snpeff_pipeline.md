### 📌 Universal End-to-End NGS Pipeline: FASTQ to Annotated VCF (SnpEff)

[![View Script](https://img.shields.io/badge/VIEW_SCRIPT-UNIVERSAL__SNPEFF__PIPELINE.PY-blue?style=for-the-badge&logo=python)](#)
[![Status](https://img.shields.io/badge/STATUS-EXECUTABLE-brightgreen?style=for-the-badge)](#)

<details>
<summary>🔍 <b>Click to View Script Details & Usage Header</b></summary>

```python
# ==============================================================================
# SCRIPT NAME : universal_snpeff_pipeline.py
# AUTHOR      : K. R. Janahan
# DESCRIPTION : Universal end-to-end NGS pipeline from raw FASTQ reads to 
#               functionally annotated and indexed VCF using BWA, Samtools, 
#               BCFtools, and an on-the-fly local SnpEff database setup.
# DEPENDENCIES: BWA, Samtools, BCFtools, SnpEff, Perl/Python3
# USAGE       :
conda activate universal_annotation_pipeline 
./universal_snpeff_pipeline.py \
  -1 your_forward_reads.fastq.gz \
  -2 your_reverse_reads.fastq.gz \
  -f your_reference_genome.fasta \
  -a your_annotation_file.gff3 \
  -g MyOrganism_v1 \
  -t 4 \
  -o final_annotated_variants.vcf.gz

# ==============================================================================
```
</details>

<details>
<summary>💻 <b>Click to View / Copy Full Source Code</b></summary>

```python
#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
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

def run_command(cmd, description):
    """Safely executes system commands via subprocess with strict error tracking."""
    print(f"\n[INFO] Starting: {description}")
    print(f"[CMD] {cmd}")
    try:
        subprocess.run(cmd, check=True, shell=True)
        print(f"[SUCCESS] Completed: {description}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed during: {description}", file=sys.stderr)
        print(f"[ERROR] Details: {e}", file=sys.stderr)
        sys.exit(1)


def setup_local_database(genome_name, fasta_path, annotation_path, is_gtf=False):
    """Generates an isolated, on-the-fly SnpEff environment and builds the database."""
    base_dir = os.path.abspath("local_snpeff_db")
    data_dir = os.path.join(base_dir, "data", genome_name)
    config_path = os.path.join(base_dir, "snpEff.config")

    print(f"\n[INFO] Setting up custom local SnpEff database for: '{genome_name}'")
    os.makedirs(data_dir, exist_ok=True)

    # 1. Write an isolated local snpEff.config
    with open(config_path, "w") as config_file:
        config_file.write(f"data.dir = {os.path.join(base_dir, 'data')}\n")
        config_file.write(f"{genome_name}.genome : {genome_name}\n")

    # 2. Link reference genome into data structure
    shutil.copyfile(fasta_path, os.path.join(data_dir, "sequences.fa"))

    # 3. Link structural annotations based on format type
    ann_filename = "genes.gtf" if is_gtf else "genes.gff"
    shutil.copyfile(annotation_path, os.path.join(data_dir, ann_filename))

    # 4. Build the SnpEff database
    fmt_flag = "-gtf22" if is_gtf else "-gff3"
    build_cmd = (
        f'snpEff build {fmt_flag} -c "{config_path}" -noCheckCds -noCheckProtein {genome_name}'
    )
    run_command(build_cmd, f"Building SnpEff DB for {genome_name}")

    return config_path


def main():
    parser = argparse.ArgumentParser(
        description="Universal End-to-End NGS Pipeline: From Raw FASTQ to Annotated VCF"
    )
    # Core Raw Inputs
    parser.add_argument("-1", "--fastq1", required=True, help="Forward raw reads (R1.fastq.gz)")
    parser.add_argument("-2", "--fastq2", required=True, help="Reverse raw reads (R2.fastq.gz)")
    parser.add_argument("-f", "--fasta", required=True, help="Reference Genome FASTA file")
    parser.add_argument("-a", "--annotation", required=True, help="Structural Genome Annotation (GFF3 or GTF)")
    parser.add_argument("-g", "--genome", required=True, help="Unique identifier name for this genome build")
    parser.add_argument("-o", "--output", required=True, help="Output Annotated VCF file path (ends in .vcf.gz)")
    parser.add_argument("--gtf", action="store_true", help="Flag: Specify if the custom annotation is GTF format (Default is GFF3)")
    parser.add_argument("-t", "--threads", default="4", help="Number of CPU threads to allocate (Default: 4)")

    args = parser.parse_args()

    # Define intermediate filenames securely
    bam_output = "sorted_alignments.bam"
    raw_vcf = "raw_variants.vcf"

    print("=========================================================")
    print("  LAUNCHING UNIVERSAL FASTQ-TO-ANNOTATED-VCF PIPELINE   ")
    print("=========================================================")

    # Step 1: Reference Genome Indexing
    if not os.path.exists(f"{args.fasta}.bwt"):
        run_command(f'bwa index "{args.fasta}"', "BWA Genome Indexing")
    if not os.path.exists(f"{args.fasta}.fai"):
        run_command(f'samtools faidx "{args.fasta}"', "Samtools Genome Indexing")

    # Step 2: Read Alignment
    alignment_pipe = (
        f'bwa mem -t {args.threads} "{args.fasta}" "{args.fastq1}" "{args.fastq2}" | '
        f'samtools view -@ {args.threads} -bS - | '
        f'samtools sort -@ {args.threads} -o "{bam_output}"'
    )
    run_command(alignment_pipe, "Paired-End Read Alignment and Coordinate Sorting")

    # Step 3: Index the BAM file
    run_command(f'samtools index "{bam_output}"', "Indexing Sorted BAM File")

    # Step 4: Variant Calling via Organism-Agnostic Bcftools
    variant_calling_cmd = (
        f'bcftools mpileup --threads {args.threads} -f "{args.fasta}" "{bam_output}" | '
        f'bcftools call --threads {args.threads} -mv -Ov -o "{raw_vcf}"'
    )
    run_command(variant_calling_cmd, "Variant Calling via bcftools")

    # Step 5: Setup Local Isolated SnpEff Database Env
    config_flag = setup_local_database(
        genome_name=args.genome,
        fasta_path=args.fasta,
        annotation_path=args.annotation,
        is_gtf=args.gtf
    )

    # Step 6: Functional Annotation & Final Packaging
    snpeff_core_cmd = f'snpEff -c "{config_flag}" {args.genome} "{raw_vcf}"'
    final_pipeline_cmd = f'{snpeff_core_cmd} | bcftools view --threads {args.threads} -Oz -o "{args.output}" && bcftools index "{args.output}"'
    run_command(final_pipeline_cmd, "SnpEff Variant Annotation and VCF Indexing")

    # Step 7: Custom Post-Processing Summaries (Dynamic f-strings properly closed!)
    print("\n=========================================================")
    print("             GENERATING POST-RUN REPORTS                 ")
    print("=========================================================")
    
    # 1. Total Line Count Output
    count_cmd = f'bcftools view -H "{args.output}" | wc -l'
    run_command(count_cmd, "Counting Final Total Variants")

    # 2. Extract Fields to annotated_variants.tsv Spreadsheet layout
    tsv_cmd = f'bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/ANN\\n" "{args.output}" > annotated_variants.tsv'
    run_command(tsv_cmd, "Exporting Formatted Variants to annotated_variants.tsv")

    # Clean Up Step
    print("\n[INFO] Cleaning up bulky intermediate alignment workflow files...")
    for template_file in [bam_output, f"{bam_output}.bai", raw_vcf]:
        if os.path.exists(template_file):
            os.remove(template_file)

    print(f"\n[SUCCESS] Pipeline completed flawlessly!")
    print(f"[OUTPUT 1] Compressed VCF Archive: {args.output}")
    print(f"[OUTPUT 2] Extracted TSV Table:     annotated_variants.tsv")
    print("=========================================================")


if __name__ == "__main__":
    main()

```
</details>
