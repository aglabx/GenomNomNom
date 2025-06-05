# 🧬 GenomNomNom v0.9

**GenomNomNom** – A playful bioinformatics tool for comprehensive genome analysis! 🍽️  
It munches through genomes, analyzes coding sequences (exons), counts all 64 codons, provides detailed statistics, and can automatically download reference genomes from NCBI.

[![Version](https://img.shields.io/badge/version-0.9-blue.svg)](https://github.com/yourusername/genomnomnom)
[![Python](https://img.shields.io/badge/python-3.8+-brightgreen.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## ✨ Features

### 🧬 Core Analysis
- ✅ **Real genome parsing** with BioPython (FASTA format)
- ✅ **Real annotation parsing** (GFF3 format) 
- ✅ **Comprehensive codon analysis** - all 64 codons with frequencies and percentages
- ✅ **Classic codon usage table** - textbook-style square table with color visualization
- ✅ **Exon (CDS) statistics** - coding sequence analysis with gene names and products
- ✅ **Start/Stop codon analysis** - detailed usage patterns
- ✅ **Genome statistics** - length, GC content, contigs

### 🌐 NCBI Integration
- ✅ **NCBI genome search** - find reference genomes by species name
- ✅ **Interactive selection** - choose from available assemblies with metadata
- ✅ **Automatic download** - both FASTA and GFF files
- ✅ **Assembly information parsing** - level, size, and quality metrics

### 📊 Output & Visualization  
- ✅ **Beautiful console output** with emojis and color coding
- ✅ **Detailed codon usage reports** with amino acid breakdowns
- ✅ **Color-coded codon table** - visual representation of usage patterns
- ✅ **CSV export** functionality for further analysis
- ✅ **Gene information** - longest/shortest exons with gene names and products

### 🛠️ User Experience
- ✅ **Command-line interface** with comprehensive help
- ✅ **Progress indicators** for long-running operations
- ✅ **Verbose mode** for detailed debugging
- ✅ **Test data included** for quick start

---

## 📦 Installation

### Requirements
- Python 3.8 or higher
- BioPython for genome parsing
- Internet connection for NCBI downloads

### Quick Install
```bash
git clone https://github.com/yourusername/genomnomnom.git
cd GenomNomNom
pip install -r requirements.txt
```

### Dependencies
```bash
pip install biopython
```

---

## 🧪 Usage

### Mode 1: Analyze Local Genome Files
```bash
# Basic analysis
python genomnomnom.py --genome genome.fasta --annotation annotation.gff

# With detailed codon table
python genomnomnom.py --genome genome.fasta --annotation annotation.gff --detailed-codons

# Save results to CSV
python genomnomnom.py --genome genome.fasta --annotation annotation.gff --output results.csv
```

### Mode 2: Search and Download from NCBI
```bash
# Interactive NCBI search
python genomnomnom.py --species "Escherichia coli" --email your.email@example.com

# With detailed codon analysis
python genomnomnom.py --species "Mycoplasma genitalium" --email your.email@example.com --detailed-codons

# Save NCBI analysis to CSV
python genomnomnom.py --species "Bacillus subtilis" --email your.email@example.com --output bacillus_results.csv
```

### Quick Start with Sample Data
```bash
# Test with included E. coli sample (5kb)
python genomnomnom.py --genome test_data/ecoli_5kb.fasta --annotation test_data/annotation_sample.gff --detailed-codons
```

---

## 📊 Example Output

### Basic Analysis Report
```
🧬 GenomNomNom Analysis Report
=====================================

📊 Genome Statistics:
- Total length: 7,887 bp
- GC content: 47.8%
- Number of contigs: 1

🚀 Start Codon Usage:
- ATG: 2 (100.0%)
- GTG: 0 (0.0%)
- TTG: 0 (0.0%)

🛑 Stop Codon Usage:
- TAA: 0 (0.0%)
- TAG: 0 (0.0%)
- TGA: 1 (100.0%)

📏 Exon Statistics:
- Total Exons: 4
- Mean length: 1187.2 bp
- Median length: 1110.0 bp
- Longest Exon: 2,463 bp
  📍 Gene: thrA
  🔬 Product: bifunctional aspartokinase I/homoserine dehydrogenase I
- Shortest Exon: 66 bp
  📍 Gene: thrL
  🔬 Product: thr operon leader peptide
```

### Detailed Codon Analysis
```
🧬 CODON USAGE BY AMINO ACID (sorted by frequency)
────────────────────────────────────────────────────────────────

🥇 🔴 Arg (R) — 8.72% of all codons
   ┌──────────────────────────────────────────────────────────────┐
   │ AGA:     38 ( 2.40% total,  27.5% of R) [█████░░░░░░░░░░░] 
   │ AGG:     23 ( 1.45% total,  16.7% of R) [███░░░░░░░░░░░░░] 
   │ CGT:     22 ( 1.39% total,  15.9% of R) [███░░░░░░░░░░░░░] 
   └──────────────────────────────────────────────────────────────┘
```

### Classic Codon Table
```
🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬
📚 CLASSIC CODON USAGE TABLE (Textbook Style)
🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬

       ┃       T        ┃       C        ┃       A        ┃       G        ┃
  1st  ┃  2nd position  ┃  2nd position  ┃  2nd position  ┃  2nd position  ┃
 pos.  ┃  T  C  A  G  ┃  T  C  A  G  ┃  T  C  A  G  ┃  T  C  A  G  ┃
━━━━━━━╋━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━╋
   T   ┃  1.3  1.0  1.3  1.8 ┃  0.6  1.5  0.7  0.9 ┃  1.3  1.4  1.1  0.6 ┃  0.8  1.3  2.4  2.7 ┃
       ┃   F   F   L   L ┃   S   S   S   S ┃   Y   Y STP STP ┃   C   C STP   W ┃

🎨 COLOR LEGEND:
• 🟢 High usage (≥2.5%)  • 🟡 Medium usage (1.5-2.4%)  • 🔵 Low usage (0.5-1.4%)
• ⚪ Very low (<0.5%)   • ⚫ Unused (0%)           • 🔴 Stop codons
```

---

## 🔧 Command Line Options

```bash
python genomnomnom.py [OPTIONS]

Required (choose one):
  --genome PATH         Path to genome FASTA file (local mode)
  --species "NAME"      Species name for NCBI search (NCBI mode)

Optional:
  --annotation PATH     Path to GFF annotation file (required for local mode)
  --email EMAIL         Email for NCBI API (required for NCBI mode)
  --output PATH         Save results to CSV file
  --detailed-codons     Show comprehensive codon usage analysis
  --verbose             Enable verbose output for debugging
  --help               Show help message and exit
```

---

## 🧬 NCBI Integration Examples

### Bacterial Genomes
```bash
# E. coli reference genome
python genomnomnom.py --species "Escherichia coli" --email your@email.com

# Minimal bacterial genome
python genomnomnom.py --species "Mycoplasma genitalium" --email your@email.com

# Model organism
python genomnomnom.py --species "Bacillus subtilis" --email your@email.com
```

### Eukaryotic Genomes
```bash
# Model yeast
python genomnomnom.py --species "Saccharomyces cerevisiae" --email your@email.com

# Human genome (warning: large download!)
python genomnomnom.py --species "Homo sapiens" --email your@email.com
```

### How NCBI Integration Works
1. 🔍 **Search**: Queries NCBI Assembly database for species
2. 📋 **Display**: Shows available assemblies with metadata (size, level, date)
3. 🎯 **Select**: Interactive choice of preferred assembly
4. 📥 **Download**: Automatic download of FASTA and GFF files
5. 🔬 **Analyze**: Full genomic analysis on downloaded data

---

## 📁 Project Structure

```
GenomNomNom/
├── 📄 genomnomnom.py          # Main analysis tool
├── 📄 requirements.txt        # Python dependencies  
├── 📄 Makefile               # Build automation
├── 📄 README.md              # This file
├── 📄 LICENSE                # MIT license
├── 📁 test_data/             # Sample genome data
│   ├── ecoli_5kb.fasta       # E. coli sample genome (5kb)
│   └── annotation_sample.gff # Corresponding annotations
├── 📁 downloads/             # NCBI downloads (auto-created)
└── 📁 drafts/                # Development notes
```

---

## 🚀 Coming in v1.0

- 🧪 **Comprehensive test suite** with pytest
- 📊 **Additional statistics** and analysis methods
- 🎨 **HTML report generation** with interactive plots
- 🔧 **Configuration file support** for analysis parameters
- 📈 **Batch processing** for multiple genomes
- 🌐 **Web interface** for browser-based analysis

---

## 💡 Tips & Best Practices

### For Local Analysis
- Ensure your GFF file matches your FASTA file (same sequence IDs)
- Use `--verbose` for troubleshooting parsing issues
- Export to CSV for integration with other bioinformatics tools

### For NCBI Downloads  
- Use a valid email address (required by NCBI Entrez)
- Choose "Complete Genome" assemblies when available
- Be aware that large genomes (e.g., human) may take time to download

### Interpreting Results
- High codon usage bias suggests strong selection pressure
- Unusual start/stop codon patterns may indicate annotation issues
- Use the detailed codon table to identify preferred codons for each amino acid

---

## 🐛 Troubleshooting

### Common Issues
- **"No CDS features found"**: Check that your GFF file contains CDS annotations
- **"Failed to download"**: Verify internet connection and email address
- **"Sequence not found"**: Ensure FASTA and GFF files match
- **Import errors**: Install BioPython with `pip install biopython`

### Getting Help
- Use `--verbose` for detailed debugging information
- Check that input files are in the correct format
- Ensure all dependencies are installed

---

## 🍽️ Why "GenomNomNom"?

Because it has an insatiable appetite for genomes! It devours FASTA files, munches through annotations, and always comes back hungry for more genomic data. 🧬🍴

---

## 🤝 Contributing

We welcome contributions! Whether it's bug reports, feature requests, or code improvements:

1. Fork the repository
2. Create a feature branch
3. Make your changes  
4. Add tests (for v1.0+)
5. Submit a pull request

---

## 📜 License

MIT License - see [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- **BioPython** for robust genomic data parsing
- **NCBI** for providing free access to genomic databases  
- **The bioinformatics community** for inspiration and feedback

---

*Made with 💙 for the genomics community*
