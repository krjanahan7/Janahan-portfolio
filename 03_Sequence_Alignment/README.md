# 🧭 03_Sequence Alignment & Mapping

![Tools](https://img.shields.io/badge/Tools-BWA%20%7C%20Bowtie2%20%7C%20HISAT2%20%7C%20SAMtools-blue)
![Language](https://img.shields.io/badge/Language-Bash-orange)
![Stage](https://img.shields.io/badge/Stage-Alignment-yellow)

---

## 📌 Objective

This section documents how I map cleaned reads to a reference genome and convert, sort, and index alignments for downstream variant calling or quantification.

---

## 🧠 Key Topics Covered

- Reference genome indexing
- DNA alignment with BWA-MEM
- Spliced RNA alignment with HISAT2
- SAM → BAM conversion, sorting, indexing
- Duplicate marking
- Alignment quality statistics

---

## 💻 Commands Used

Index the reference:

    bwa index reference.fasta
    samtools faidx reference.fasta

Align DNA reads:

    bwa mem -t 4 -R "@RG\tID:S1\tSM:sample1\tPL:ILLUMINA" \
      reference.fasta sample_R1.fq.gz sample_R2.fq.gz > sample.sam

Convert, sort, index:

    samtools view -bS sample.sam > sample.bam
    samtools sort -@ 4 -o sample.sorted.bam sample.bam
    samtools index sample.sorted.bam

Alignment statistics:

    samtools flagstat sample.sorted.bam
    samtools depth -a sample.sorted.bam | awk '{s+=$3} END {print s/NR}'

Spliced alignment (RNA-Seq):

    hisat2-build reference.fasta genome_index
    hisat2 -p 4 -x genome_index -1 R1.fq.gz -2 R2.fq.gz -S rna.sam

---

## ⚙️ Parameter Explanations

| Parameter | Meaning |
|---|---|
| `-t 4` / `-p 4` | Number of CPU threads |
| `-R "@RG..."` | Read group tag required by GATK |
| `-bS` | Input SAM, output BAM |
| `sort` | Orders alignments by coordinate |
| `index` | Creates `.bai` for fast random access |
| `flagstat` | Summarises mapped, paired, and duplicate reads |

---

## 📊 Results & Interpretation

- Overall mapping rate above ~90% indicates a suitable reference and clean reads.
- Properly paired percentage reflects library integrity.
- Mean depth guides confidence in variant calls.

---

## ⭐ Key Takeaway

Alignment quality sets the ceiling for every downstream result — a poorly mapped BAM cannot be rescued by better variant callers.

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
