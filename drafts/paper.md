# GenomNomNom: A Playful Toolkit for Rapid Exploratory Gene‑Structure Analysis Across Diverse Genomes

**Authors**: Special Topics in Genomics course participants\*, ITMO University, St Petersburg, Russia
*All authors contributed equally to this draft.*

---

## Abstract

GenomNomNom is an open‑source Python toolkit that “munches” through genomes to deliver fast, interpretable statistics on gene architecture.  In its initial release the tool ingests a genome FASTA file and a matching annotation and outputs concise summaries of start/stop codon usage and open reading frame (ORF) distribution.  Designed for education and rapid hypothesis generation, GenomNomNom lowers the entry barrier to comparative gene‑structure analysis and lays the groundwork for deeper, HMM‑free gene prediction.  Here we describe the motivation, core design, and planned extensions, and outline a set of teaching case studies that will accompany future versions.

**Keywords**: ORF, codon usage, gene prediction, comparative genomics, teaching resource, Python toolkit

---

## 1 Introduction

The accurate identification of protein‑coding genes remains central to genome annotation.  Mature pipelines—e.g. Glimmer 3, Prodigal, GeneMark—rely on Hidden Markov Models (HMMs) and species‑specific training, creating a steep learning curve for students and complicating exploratory analyses across multiple taxa.  Meanwhile, biologically illuminating questions can often be asked (and partly answered) with simpler statistics: How balanced are start codons?  Do stop‑codon preferences shift in reduced genomes?  How many frames remain unutilised?  GenomNomNom was conceived as a lightweight answer to these questions, suitable for classroom use yet extensible to research‑grade workflows.

## 2 Design and Implementation

### 2.1 Architecture overview

* **Input**: Genome in FASTA format; gene annotation in GFF3 (or compatible) format.
* **Core engine**: Pandas‑based parser that maps coding regions, derives per‑codon counts, and summarises ORFs.
* **Output**: Markdown/CSV table with counts and proportions of each canonical start (ATG/GTG/TTG) and stop (TAA/TAG/TGA) codon; ORF length distribution metrics.

### 2.2 Extensibility roadmap

1. Built‑in ORF scanner for unannotated genomes (threshold configurable).
2. Comparative module aligning codon statistics across species.
3. Plug‑in API for external predictors (e.g. Glimmer, GeneMark) enabling head‑to‑head benchmarking.
4. Interactive visualisations (Plotly/Altair) and Jupyter notebook templates.

## 3 Materials & Methods

### 3.1 Data sets (planned)

* *Escherichia coli* K‑12 MG1655 (RefSeq NC\_000913.3)
* *Mycoplasma genitalium* G37 (RefSeq NC\_000908.2)
* Additional draft genomes representing AT‑ and GC‑rich extremes (TBD)

### 3.2 Statistical analysis

Counts are normalised to coding sequence length (codons per kilobase).  Confidence intervals for codon proportions will be computed using Wilson score intervals.  Cross‑species differences will be assessed via χ² tests with Bonferroni correction.

## 4 Preliminary Results *(placeholders)*

| Genome          | ATG (%) | GTG (%) | TTG (%) | TAA (%) | TAG (%) | TGA (%) | Mean ORF length (bp) |
| --------------- | ------- | ------- | ------- | ------- | ------- | ------- | -------------------- |
| *E. coli*       | TBD     | TBD     | TBD     | TBD     | TBD     | TBD     | TBD                  |
| *M. genitalium* | TBD     | TBD     | TBD     | TBD     | TBD     | TBD     | TBD                  |

> Detailed results will be added after the first analysis run.

## 5 Discussion and Future Work

Our initial statistics will provide a baseline for exploring codon‑usage evolution and testing hypotheses such as the enrichment of TGA‑>Trp recoding in obligate parasites.  Future iterations will integrate lightweight machine‑learning classifiers for de‑novo gene prediction and a benchmark suite against established HMM‑based tools.

## 6 Availability and Requirements

* **Project home page**: [https://github.com/yourusername/GenomNomNom](https://github.com/aglabx/GenomNomNom)
* **Operating system**: Platform‑independent
* **Programming language**: Python ≥3.10
* **License**: MIT
* **Dependencies**: pandas, biopython, click (full list in `requirements.txt`)

## 7 Acknowledgements

We thank the instructors of the “Special Topics in Genomics” course at ITMO University for guidance and inspiration.

## 8 References *(to be completed)*

1. Delcher AL *et al.* Fast algorithms for large‑scale genome alignment and comparison. *Nucleic Acids Res.* 2002.
2. Hyatt D *et al.* Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics.* 2010.
3. Borodovsky M, Lomsadze A. GeneMark: Web software for gene finding in prokaryotes, eukaryotes and viruses. *Nucleic Acids Res.* 2011.
