### 📌 Universal Transcriptomics Pipeline: FASTQ to Gene-Level Count Matrix

[![View Script](https://img.shields.io/badge/VIEW_SCRIPT-UNIVERSAL__TRANSCRIPTOMICS__PIPELINE.PY-blue?style=for-the-badge&logo=python)](#)
[![Status](https://img.shields.io/badge/STATUS-EXECUTABLE-brightgreen?style=for-the-badge)](#)

<details>
<summary>🔍 <b>Click to View Script Details & Usage Header</b></summary>

```python
# ==============================================================================
# SCRIPT NAME : universal_transcriptomics_pipeline.py
# AUTHOR      : K. R. Janahan
# DESCRIPTION : Universal RNA-seq quantification pipeline from raw FASTQ reads 
#               to gene-level count matrices using HISAT2, Samtools, and Subread featureCounts.
# DEPENDENCIES: HISAT2, Samtools, featureCounts (Subread), Python3
# USAGE       :
conda activate universal_transcriptomics_pipeline
chmod +x universal_transcriptomics_pipeline.py
./universal_transcriptomics_pipeline.py \
  --ref_fna <YOUR_REFERENCE.fna> \
  --ref_gff <YOUR_ANNOTATION.gff> \
  --read1 <NEW_SAMPLE_1.fastq.gz> \
  --read2 <NEW_SAMPLE_2.fastq.gz> \
  --prefix <NEW_SAMPLE_NAME> \
  --threads 16

# ==============================================================================
```
</details>

<details>
<summary>💻 <b>Click to View / Copy Full Source Code</b></summary>

```python

#!/usr/bin/env python3
"""Universal Transcriptomics Pipeline — FASTQ to gene-level count matrix.

Usage:
    ./universal_transcriptomics_pipeline.py \
        --ref_fna REF.fna[.gz] \
        --ref_gff REF.gff[.gz] \
        --read1 READS_1.fastq[.gz] \
        --read2 READS_2.fastq[.gz] \
        --prefix SAMPLE_PREFIX \
        --threads N
"""

import argparse
import gzip
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

REQUIRED_TOOLS = ["hisat2-build", "hisat2", "samtools", "featureCounts"]


def log(msg):
    print(f"[PIPELINE] {msg}", flush=True)


def run(cmd):
    log("RUN: " + " ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        log(f"ERROR: '{cmd[0]}' failed (exit {result.returncode})")
        sys.exit(result.returncode)


def check_tools():
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        log("ERROR: missing tools: " + ", ".join(missing))
        log("Install with: conda install -c bioconda " + " ".join(missing))
        sys.exit(1)


def check_inputs(paths):
    missing = [p for p in paths if not os.path.isfile(p)]
    if missing:
        log("ERROR: input file(s) not found:")
        for m in missing:
            log(f"   - {m}")
        sys.exit(1)


def decompress(path):
    """Decompress .gz inputs to a sibling file, leaving the original intact."""
    if not path.endswith(".gz"):
        return path, False
    out_path = path[:-3]
    if not os.path.isfile(out_path):
        log(f"Decompressing {path} -> {out_path}")
        with gzip.open(path, "rb") as f_in, open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return out_path, True


def main():
    p = argparse.ArgumentParser(description="Universal RNA-seq quantification pipeline")
    p.add_argument("--ref_fna", required=True)
    p.add_argument("--ref_gff", required=True)
    p.add_argument("--read1", required=True)
    p.add_argument("--read2", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--keep-temp", action="store_true",
                    help="Keep decompressed inputs, SAM, and unsorted BAM")
    args = p.parse_args()

    # 1. Pre-flight: sanity checks
    check_tools()
    check_inputs([args.ref_fna, args.ref_gff, args.read1, args.read2])

    ref_fna, t1 = decompress(args.ref_fna)
    ref_gff, t2 = decompress(args.ref_gff)
    read1, t3 = decompress(args.read1)
    read2, t4 = decompress(args.read2)
    temp_files = [f for f, made in
                   [(ref_fna, t1), (ref_gff, t2), (read1, t3), (read2, t4)] if made]

    # 2. Indexing
    index_dir = f"{args.prefix}_hisat2_index"
    os.makedirs(index_dir, exist_ok=True)
    index_prefix = os.path.join(index_dir, args.prefix)
    run(["hisat2-build", "-p", str(args.threads), ref_fna, index_prefix])

    # 3. Alignment
    sam_file = f"{args.prefix}.sam"
    summary_file = f"{args.prefix}_alignment_summary.txt"
    run(["hisat2", "-p", str(args.threads),
         "-x", index_prefix,
         "-1", read1, "-2", read2,
         "-S", sam_file,
         "--summary-file", summary_file])

    # 4. SAM -> BAM -> sorted -> indexed
    bam_file = f"{args.prefix}.bam"
    sorted_bam = f"{args.prefix}_sorted.bam"
    run(["samtools", "view", "-@", str(args.threads), "-b", sam_file, "-o", bam_file])
    run(["samtools", "sort", "-@", str(args.threads), bam_file, "-o", sorted_bam])
    run(["samtools", "index", "-@", str(args.threads), sorted_bam])

    # 5. Quantification
    counts_raw = f"{args.prefix}_counts.txt"
    run(["featureCounts", "-T", str(args.threads), "-p",
         "-t", "gene", "-g", "locus_tag",
         "-a", ref_gff, "-o", counts_raw, sorted_bam])

    # 6. Simplified Gene ID + Count TSV
    counts_tsv = f"{args.prefix}_counts.tsv"
    with open(counts_raw) as f_in, open(counts_tsv, "w") as f_out:
        lines = f_in.readlines()
        f_out.write("Geneid\tCount\n")
        for line in lines[2:]:
            cols = line.rstrip("\n").split("\t")
            f_out.write(f"{cols[0]}\t{cols[-1]}\n")

    # Cleanup intermediates (originals are never touched)
    if not args.keep_temp:
        for f in temp_files + [sam_file, bam_file]:
            if os.path.isfile(f):
                os.remove(f)

    log("Done.")
    log(f"  Sorted BAM    : {sorted_bam}")
    log(f"  Raw counts    : {counts_raw}")
    log(f"  Count matrix  : {counts_tsv}")


if __name__ == "__main__":
    main()

```
</details>
