# GenomNomNom: A Playful Toolkit for Genome‑Scale Codon and ORF Exploration

**Authors**: Special Topics in Genomics course participants\*, ITMO University, St Petersburg, Russia
*All authors contributed equally to this draft.*

---

## Abstract

**GenomNomNom** is an open‑source Python toolkit that “munches” through public or local genomes to deliver instant, interpretable statistics on codon usage and gene architecture.  Given a species name—or a local FASTA + GFF pair—the program automatically fetches the reference assembly from NCBI, parses the annotation, and produces a comprehensive report covering every codon (64/64), start/stop preferences, open reading‑frame (ORF) length distribution, and basic genome descriptors.  Designed for education and rapid hypothesis generation, GenomNomNom lowers the entry barrier to comparative codon‑bias analysis while providing hooks for deeper, HMM‑free gene prediction.  Here we describe the motivation, implementation, and planned extensions, and outline benchmark scenarios that will accompany future releases.

**Keywords**: codon usage, ORF, gene prediction, comparative genomics, teaching resource, Python toolkit, NCBI API

---

## 1 Introduction

Codon‑usage bias and start/stop codon choice vary across taxa and can illuminate evolutionary forces ranging from mutational pressure to translational selection.  Established pipelines (e.g. Glimmer 3, Prodigal, GeneMark) rely on Hidden Markov Models (HMMs) and organism‑specific training, creating a steep learning curve for newcomers.  Simple summary statistics, however, already answer many biologically relevant questions: Which stop codon dominates in AT‑rich parasites?  Does TGA encode tryptophan in close relatives of *Mycoplasma*?  How skewed are codon preferences in nematodes versus streptococci?  GenomNomNom was conceived as an ultra‑lightweight answer.  A single command fetches a reference genome, counts every codon, and returns a human‑readable report ready for classroom discussion or exploratory research.

## 2 Design and Implementation

\### 2.1 Command‑line interface

```bash
python genomnomnom.py \
    --species "Caenorhabditis elegans" \
    --email user@example.com \
    --detailed-codons
```

*If files are already local, supply `--genome` and `--annotation` instead of `--species`.*

\### 2.2 Workflow

1. **NCBI fetcher** – Queries the Assembly database, lists available references, downloads FASTA and GFF.
2. **Genome parser** – Streams sequences, reports size, GC content, contig count.
3. **Annotation parser** – Extracts CDS features, merges split locations, records gene IDs/products.
4. **Codon counter** – Tallies all 64 codons across CDS; separately counts start and stop codons.
5. **ORF analyzer** – Computes length distribution, longest/shortest ORFs, and frame statistics.
6. **Report generator** – Renders markdown/plain‑text tables, optional CSV/JSON output for downstream analysis.

\### 2.3 Extensibility roadmap

* Built‑in de‑novo ORF scanner for unannotated genomes.
* Comparative module aligning codon statistics across arbitrary species sets (e.g. Gram‑positive vs. Gram‑negative).
* API wrappers for Glimmer/Prodigal to enable head‑to‑head benchmarking.
* Interactive visualisation notebooks (Plotly, Altair).
* Snakemake workflow for large‑scale automated surveys.

## 3 Materials & Methods (planned)

\### 3.1 Data sets
Initial demonstrations will use:

* *Caenorhabditis elegans* (Eukaryota; reference WBcel235)
* *Escherichia coli* K‑12 MG1655 (Proteobacteria; RefSeq NC\_000913.3)
* *Mycoplasma genitalium* G37 (Tenericutes; RefSeq NC\_000908.2)
* A GC‑rich Gram‑positive representative (*Micrococcus luteus*) and a GC‑poor Gram‑negative draft (TBD).

\### 3.2 Codon‑usage statistics
Counts are normalised per 1 000 codons; amino‑acid–specific frequency tables follow Sharp & Li’s codon‑adaptation framework.  Effective number of codons (Nc) and relative synonymous codon usage (RSCU) will be calculated; differences between groups assessed by principal‑component analysis (PCA) and MANOVA.

## 4 Preliminary Example Output

> **Command:** `genomnomnom.py --species "Caenorhabditis elegans" --detailed-codons`
>
> **Summary (excerpt):**
> • Genome length = 100 Mb, GC = 35.4 %
> • Genes = 44 795; CDS = 204 522
> • Start codons: ATG 85.9 %, GTG 9.8 %, TTG 4.3 %
> • Stop codons: TAA 43.8 %, TAG 25.4 %, TGA 30.8 %
> • Total codons analysed: 6 008 101 (all 64 in use)
> • Most frequent amino acid: Leu 8.74 %; least: Trp 1.28 %
> • Codon‑usage bias ratio (max/min): 13.2 : 1

A complete codon table is provided as supplementary CSV.

## 5 Discussion and Future Work

The automated species‑name workflow fills a gap between heavy annotation suites and manual coding.  By simplifying large‑scale codon surveys—e.g., contrasting Gram‑positive and Gram‑negative taxa—GenomNomNom enables rapid hypothesis testing in evolutionary genomics and microbiology pedagogy.  Planned de‑novo ORF detection will further position the tool as a lightweight alternative to HMM‑based predictors, while comparative visual dashboards will appeal to bench scientists exploring newly sequenced isolates.

## 6 Availability and Requirements

* **Source code**: [https://github.com/yourusername/GenomNomNom](https://github.com/yourusername/GenomNomNom)
* **OS**: Platform‑independent
* **Language**: Python ≥ 3.10
* **License**: MIT
* **Dependencies**: pandas, biopython, click, rich, requests

## 7 Acknowledgements

We thank the instructors of the “Special Topics in Genomics” course at ITMO University for guidance, and the NCBI help‑desk for API support.

## 8 References *(to be completed)*

1. Sharp PM, Li WH (1987) The codon Adaptation Index—A measure of directional synonymous codon usage bias. *Nucleic Acids Res.*
2. Delcher AL *et al.* Fast algorithms for large‑scale genome alignment and comparison. *Nucleic Acids Res.* 2002.
3. Hyatt D *et al.* Prodigal: prokaryotic gene recognition and translation‑initiation site identification. *BMC Bioinformatics* 2010.
4. Borodovsky M, Lomsadze A (2011) GeneMark: gene prediction with automated training for prokaryotic genomes. *Nucleic Acids Res.*
