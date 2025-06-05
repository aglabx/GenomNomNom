# 🧬 GenomNomNom v0.9

**GenomNomNom** – A comprehensive genomic analysis tool that combines real genome parsing, codon usage analysis, exon statistics, and NCBI integration. Perfect for researchers, students, and bioinformatics enthusiasts who want to "munch through genomes" with st---

## 🍽️ Why "GenomNomNom"?

Because it hungrily devours genomes and spits out delicious insights! This tool was born from the idea that genomic analysis should be both powerful and fun. Whether you're a researcher analyzing bacterial genomes, a student learning bioinformatics, or just curious about the genetic code, GenomNomNom makes genome analysis accessible and enjoyable.

The playful name reflects our philosophy: science doesn't have to be serious all the time. Sometimes the best insights come when you're having fun exploring data! 🧬🍽️

---

## 📜 License

MIT License - feel free to use, modify, and distribute!

---

## 📞 Support

- **Issues**: Report bugs and request features on GitHub
- **Questions**: Start a discussion in the GitHub repository  
- **Email**: Contact maintainers for urgent issues

**Version**: 0.9 (Ready for serious genomic analysis!)
**Status**: Production ready, comprehensive testing planned for v1.0--

## 🚀 Features

### ✅ **Core Analysis Capabilities**
- **Real genome parsing** with BioPython (FASTA format)
- **Real annotation parsing** (GFF3 format) 
- **Comprehensive codon usage analysis** with frequency counting and bias detection
- **Classical codon table visualization** with color-coded frequency levels
- **Exon analysis** including length distribution and gene product identification
- **Start/stop codon statistics** with detailed usage patterns

### ✅ **NCBI Integration** 
- **Automated genome search** by species name
- **Interactive genome selection** from search results with metadata
- **Assembly level filtering** (Complete, Chromosome, Scaffold, Contig)
- **Automatic download** of both FASTA and GFF files
- **Robust error handling** for network issues

### ✅ **Output & Visualization**
- **Beautiful console output** with emojis and color coding
- **CSV export functionality** for further analysis
- **Detailed statistical reports** with comprehensive genome metrics
- **Classical textbook-style codon tables** with ANSI color visualization
- **Gene product information** for longest/shortest exons

### ✅ **User Experience**
- **Command-line interface** with comprehensive help
- **Verbose mode** for detailed debugging output
- **Progress indicators** for long-running operations
- **Error recovery** and informative error messages

---

## 📦 Installation

### Prerequisites
- Python 3.7+ (tested with Python 3.11)
- pip package manager

### Install from source
```bash
git clone https://github.com/yourusername/genomnomnom.git
cd genomnomnom
pip install -r requirements.txt
```

### Dependencies
- `biopython` - For genome and annotation parsing
- `requests` - For NCBI API communication  
- `argparse` - Command-line interface (built-in)`

---

## 🧪 Usage

GenomNomNom supports two main analysis modes:

### Mode 1: Analyze Local Genome Files
```bash
python genomnomnom.py --genome genome.fasta --annotation annotation.gff
```

### Mode 2: Search and Download from NCBI
```bash
python genomnomnom.py --species "Escherichia coli" --email your.email@example.com
```

### Additional Options

#### Save results to CSV:
```bash
# Local files
python genomnomnom.py --genome genome.fasta --annotation annotation.gff --output results.csv

# NCBI download
python genomnomnom.py --species "Mycoplasma genitalium" --email your.email@example.com --output myco_results.csv
```

#### Enable verbose output for debugging:
```bash
python genomnomnom.py --genome genome.fasta --annotation annotation.gff --verbose
```

#### Quick test with sample data:
```bash
# Test local file analysis
python genomnomnom.py --genome test_data/ecoli_5kb.fasta --annotation test_data/annotation_sample.gff

# Test NCBI functionality (requires internet)
python demo_ncbi.py
```

### Using the Makefile (Advanced):
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

---

## 📊 Analysis Output

GenomNomNom provides comprehensive genomic analysis including:

### 🧬 **Genome Statistics**
- Total genome length and GC content
- Number of contigs/chromosomes
- Assembly level and metadata (for NCBI genomes)

### 🔤 **Codon Usage Analysis**
- Complete codon frequency counting from CDS sequences
- Start codon usage patterns (ATG, GTG, TTG, etc.)
- Stop codon distribution (TAA, TAG, TGA)
- Codon bias detection and statistics

### 📚 **Classical Codon Table**
A beautiful textbook-style 4x4x4 codon table with color-coded frequency visualization:
- 🔴 **Red**: High frequency codons (>75th percentile)
- 🟡 **Yellow**: Medium frequency codons (25th-75th percentile)  
- 🔵 **Blue**: Low frequency codons (<25th percentile)
- ⚫ **Gray**: Unused codons

```
📚 CLASSIC CODON USAGE TABLE (Textbook Style)

     T        C        A        G
T  TTT(12)  TCT(8)   TAT(15)  TGT(3)   T
   TTC(18)  TCC(12)  TAC(22)  TGC(7)   C  
   TTA(5)   TCA(6)   TAA(0)*  TGA(0)*  A
   TTG(14)  TCG(9)   TAG(0)*  TGG(11)  G

(Numbers show usage frequency, * indicates stop codons)
```

### 🧩 **Exon Analysis**
- Total number of coding sequence exons
- Length distribution statistics (min, max, mean, median)
- Gene products for longest and shortest exons
- Exon count per gene analysis

### 📈 **Export Options**
- Console output with beautiful formatting and emojis
- CSV export for spreadsheet analysis and further processing
- Detailed verbose logging for troubleshooting

---

## 🌐 NCBI Integration

GenomNomNom seamlessly integrates with NCBI's databases for automated genome retrieval:

### Search Process
1. **Query NCBI Assembly database** by species name
2. **Filter results** by assembly level (Complete, Chromosome, Scaffold, Contig)
3. **Display interactive selection** with metadata:
   - Assembly accession and name
   - Genome size and assembly level
   - Submission date and organism info
4. **Automatic download** of FASTA and GFF files
5. **Run complete analysis** on downloaded data

### Example Commands
```bash
# Search for E. coli reference genomes
python genomnomnom.py --species "Escherichia coli" --email your.email@example.com

# Search for human reference genome
python genomnomnom.py --species "Homo sapiens" --email your.email@example.com

# Search for a specific bacterial strain
python genomnomnom.py --species "Mycoplasma genitalium" --email your.email@example.com
```

### Features
- **Robust error handling** for network connectivity issues
- **Assembly level filtering** (prioritizes Complete > Chromosome > Scaffold > Contig)
- **Interactive selection** from multiple matching genomes
- **Automatic file organization** in `downloads/` directory
- **Progress indicators** for download operations

---

## 🛠️ Troubleshooting

### Common Issues

**"No CDS features found"**: Your GFF file may not contain CDS (Coding Sequence) annotations. GenomNomNom requires CDS features for codon and exon analysis.

**"Network timeout"**: NCBI servers may be busy. Try again later or use `--verbose` flag to see detailed error messages.

**"Email required for NCBI"**: NCBI requires an email address for API access. Use `--email your.email@example.com`.

**"BioPython import error"**: Install dependencies with `pip install -r requirements.txt`.

### Getting Help
- Use `--help` flag for command-line options
- Use `--verbose` flag for detailed debugging output
- Check the `test_data/` directory for sample files
- Run `python demo_ncbi.py` to test NCBI functionality

---

## � Project Structure

```
GenomNomNom/
├── genomnomnom.py              # Main analysis tool with all functionality
├── requirements.txt            # Python dependencies
├── Makefile                   # Build automation and shortcuts
├── README.md                  # This documentation
├── LICENSE                    # MIT License
├── test_data/                 # Sample genome data for testing
│   ├── ecoli_5kb.fasta       # E. coli genome segment
│   └── annotation_sample.gff  # Matching GFF annotation
├── downloads/                 # NCBI downloaded genomes (auto-created)
├── drafts/                    # Development documentation
└── test_*.py                  # Test scripts (not updated for v0.9)
```

### Core Classes
- **NCBISearcher**: Handles NCBI database queries and downloads
- **GenomeParser**: FASTA file parsing and sequence analysis  
- **AnnotationParser**: GFF3 file parsing and feature extraction
- **CodonAnalyzer**: Codon usage analysis and table generation
- **ExonAnalyzer**: Coding sequence exon analysis
- **ReportGenerator**: Output formatting and CSV export

---

## 🔮 Roadmap to v1.0

### Planned Features
- **Comprehensive test suite** with unit and integration tests
- **Advanced filtering options** for NCBI search results
- **Comparative genomics** between multiple species
- **Interactive web interface** with Streamlit/Flask
- **HMM-based gene prediction** for unannotated genomes
- **Phylogenetic analysis** integration
- **Performance optimizations** for large genomes

### Quality Improvements
- **Input validation** and error handling enhancements
- **Memory optimization** for large genome processing
- **Parallel processing** support for multi-genome analysis
- **Configuration file** support for batch operations
- **Docker containerization** for easy deployment

---

## 🤝 Contributing

GenomNomNom is open source and welcomes contributions! Areas where help is needed:

- **Testing**: Writing comprehensive test suites
- **Documentation**: Improving examples and tutorials
- **Features**: Implementing advanced genomic analysis features
- **Performance**: Optimizing for large-scale genomic data
- **UI/UX**: Creating web interfaces and better visualizations

---

## 🙏 Acknowledgments

- **BioPython community** for excellent genomic data parsing tools
- **NCBI** for providing comprehensive genomic databases and APIs
- **Python ecosystem** for making bioinformatics accessible
- **Open source community** for inspiration and best practices

Special thanks to researchers and students who provided feedback during development!

---

## 🍽️ Why "GenomNomNom"?

Because it eats genomes. And it’s hungry for more.

---

## 📜 License

MIT License

```
