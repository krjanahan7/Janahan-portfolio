# 🧩 07_Genome Assembly & Annotation

![Tools](https://img.shields.io/badge/Tools-SPAdes%20%7C%20QUAST%20%7C%20Prokka%20%7C%20BUSCO-blue)
![Language](https://img.shields.io/badge/Language-Bash-orange)
![Stage](https://img.shields.io/badge/Stage-Assembly-yellow)

---

## 📌 Objective

This section documents de novo assembly of sequencing reads into contigs, evaluation of assembly quality, and structural/functional genome annotation.

---

## 🧠 Key Topics Covered

- De Bruijn graph assembly concepts and k-mer selection
- De novo assembly with SPAdes
- Assembly metrics: N50, L50, total length, GC content
- Completeness assessment with BUSCO
- Prokaryotic annotation with Prokka

---

## 💻 Commands Used

Assemble:

    spades.py -1 R1_paired.fq.gz -2 R2_paired.fq.gz \
              -o spades_output -t 4 -m 16 --careful

Evaluate:

    quast.py spades_output/contigs.fasta -o quast_report
    busco -i spades_output/contigs.fasta -l bacteria_odb10 -m genome -o busco_out

Annotate:

    prokka --outdir prokka_out --prefix sample spades_output/contigs.fasta

---

## ⚙️ Parameter Explanations

| Parameter | Meaning |
|---|---|
| `--careful` | Reduces mismatches and short indels in the assembly |
| `-m 16` | Memory limit in GB |
| `-l bacteria_odb10` | BUSCO lineage dataset for completeness scoring |
| `-m genome` | BUSCO assessment mode |
| `--prefix` | Naming prefix for annotation output files |

---

## 📊 Results & Interpretation

- Higher N50 with fewer contigs indicates a more contiguous assembly.
- BUSCO completeness above ~95% suggests most conserved genes are captured.
- Prokka outputs GFF, GBK, and FAA files listing predicted CDS, tRNA, and rRNA features.

---

## ⭐ Key Takeaway

Assembly quality is judged by contiguity *and* completeness — a long N50 means little if conserved genes are missing.

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
