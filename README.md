# 🧬 GenomNomNom

**GenomNomNom** – A playful tool for munching through genomes: finds genes, counts start/stop codons, analyzes ORFs, and compares annotation strategies.

---

## 🚀 Features (current mock version)

- ✅ Parses a genome in FASTA format (file validation)
- ✅ Parses gene annotations (GFF format validation)
- ✅ Mock codon counting functionality:
  - Start codons (ATG, GTG, TTG) with percentages
  - Stop codons (TAA, TAG, TGA) with percentages
- ✅ Mock ORF statistics (length distribution)
- ✅ Beautiful console output with emojis
- ✅ CSV export functionality
- ✅ Command-line interface with help
- ✅ Test data and automated testing

**Note:** This is currently a mock implementation that demonstrates the expected interface and output format. The actual genome parsing and analysis will be implemented in future versions.

---

## 📦 Installation

```bash
git clone https://github.com/yourusername/genomnomnom.git
cd genomnomnom
pip install -r requirements.txt
````

---

## 🧪 Usage

### Basic usage:
```bash
python genomnomnom.py --genome genome.fasta --annotation annotation.gff
```

### Save results to CSV:
```bash
python genomnomnom.py --genome genome.fasta --annotation annotation.gff --output results.csv
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
python genomnomnom.py --genome test_data/genome_sample.fasta --annotation test_data/annotation_sample.gff
```

Outputs a comprehensive analysis report including:
- Genome statistics (length, GC content, contigs)
- Start/stop codon usage with percentages  
- ORF length distribution statistics

---

## 🔮 Coming soon

* ORF detection without annotation
* HMM-free gene prediction
* Comparative genomics
* Interactive visualizations

---

## 🍽️ Why "GenomNomNom"?

Because it eats genomes. And it’s hungry for more.

---

## 📜 License

MIT License

```
