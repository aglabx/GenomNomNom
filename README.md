# 🧬 GenomNomNom

**GenomNomNom** – A playful tool for munching through genomes: finds genes, counts start/stop codons, analyzes ORFs, and can automatically download reference genomes from NCBI.

---

## 🚀 Features

- ✅ **Real genome parsing** with BioPython (FASTA format)
- ✅ **Real annotation parsing** (GFF3 format)
- ✅ **NCBI integration**: Search and download reference genomes by species name
- ✅ **Interactive genome selection** from NCBI search results
- ✅ Mock codon counting functionality (will be replaced with real analysis)
- ✅ Mock ORF statistics (will be replaced with real analysis)
- ✅ Beautiful console output with emojis
- ✅ CSV export functionality
- ✅ Command-line interface with help
- ✅ Test data and automated testing

---

## 📦 Installation

```bash
git clone https://github.com/yourusername/genomnomnom.git
cd genomnomnom
pip install -r requirements.txt
````

---

## 🧪 Usage

### Mode 1: Analyze local genome files
```bash
python genomnomnom.py --genome genome.fasta --annotation annotation.gff
```

### Mode 2: Search and download from NCBI
```bash
python genomnomnom.py --species "Escherichia coli" --email your.email@example.com
```

### Save results to CSV:
```bash
# Local files
python genomnomnom.py --genome genome.fasta --annotation annotation.gff --output results.csv

# NCBI download
python genomnomnom.py --species "Mycoplasma genitalium" --email your.email@example.com --output myco_results.csv
```

### With verbose output:
```bash
python genomnomnom.py --genome genome.fasta --annotation annotation.gff --verbose
```

### Using the Makefile:
```bash
# Install dependencies
make install

# Run demo with sample data
make demo

# Run tests
make test

# Clean temporary files
make clean
```

### Quick test with sample data:
```bash
# Test local file analysis
python genomnomnom.py --genome test_data/genome_sample.fasta --annotation test_data/annotation_sample.gff

# Test NCBI search (non-interactive test)
python test_ncbi.py
```

Outputs a comprehensive analysis report including:
- Genome statistics (length, GC content, contigs)
- Start/stop codon usage with percentages  
- ORF length distribution statistics

### NCBI Integration Examples:

The tool can automatically search NCBI for reference genomes and download them:

```bash
# Search for E. coli reference genomes
python genomnomnom.py --species "Escherichia coli" --email your.email@example.com

# Search for human reference genome
python genomnomnom.py --species "Homo sapiens" --email your.email@example.com

# Search for a specific bacterial strain
python genomnomnom.py --species "Mycoplasma genitalium" --email your.email@example.com
```

The tool will:
1. Search NCBI's Assembly database for reference genomes
2. Display available options with metadata (size, level, date)
3. Let you interactively select which genome to download
4. Download both FASTA and GFF files
5. Run the full analysis on the downloaded data

---

## 🔮 Coming soon

* **Real codon counting** based on CDS sequences
* **Real ORF detection** without annotation  
* HMM-free gene prediction
* Comparative genomics
* Interactive visualizations
* Advanced NCBI filtering options

---

## 📁 Project Structure

```
GenomNomNom/
├── genomnomnom.py          # Main analysis tool
├── test_ncbi.py            # NCBI functionality test
├── test_genomnomnom.py     # General tests
├── requirements.txt        # Python dependencies
├── Makefile               # Build automation
├── test_data/             # Sample genome data
├── downloads/             # Downloaded NCBI genomes
└── docs/                  # Documentation
```

---

## 🍽️ Why "GenomNomNom"?

Because it eats genomes. And it’s hungry for more.

---

## 📜 License

MIT License

```
