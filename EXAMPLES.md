# GenomNomNom Usage Examples

## 🚀 Quick Start

### 1. Basic Analysis
```bash
python genomnomnom.py --genome test_data/genome_sample.fasta --annotation test_data/annotation_sample.gff
```

Expected output:
```
🧬 GenomNomNom Analysis Report
=====================================

📊 Genome Statistics:
- Total length: 4,641,652 bp
- GC content: 50.8%
- Number of contigs: 1

🚀 Start Codon Usage:
- ATG: 3,847 (95.6%)
- GTG: 152 (3.8%)
- TTG: 24 (0.6%)

🛑 Stop Codon Usage:
- TAA: 1,834 (45.6%)
- TAG: 289 (7.2%)
- TGA: 1,900 (47.2%)

📏 ORF Statistics:
- Total ORFs: 4,023
- Mean length: 924.3 bp
- Median length: 678 bp
- Longest ORF: 4,932 bp
- Shortest ORF: 96 bp

🍽️ Nom nom nom! Analysis complete.
```

### 2. Save Results to CSV
```bash
python genomnomnom.py --genome test_data/genome_sample.fasta --annotation test_data/annotation_sample.gff --output my_results.csv
```

This creates a CSV file with all statistics:
```csv
metric_type,codon,count
start_codons,ATG,3847.0
start_codons,GTG,152.0
start_codons,TTG,24.0
stop_codons,TAA,1834.0
stop_codons,TAG,289.0
stop_codons,TGA,1900.0
orf_stats,total_orfs,4023.0
orf_stats,mean_length,924.3
orf_stats,median_length,678.0
orf_stats,longest_orf,4932.0
orf_stats,shortest_orf,96.0
```

### 3. Verbose Mode
```bash
python genomnomnom.py --genome test_data/genome_sample.fasta --annotation test_data/annotation_sample.gff --verbose
```

Shows detailed progress information during analysis.

## 🛠️ Using the Makefile

### Install Dependencies
```bash
make install
```

### Run Demo
```bash
make demo
```

### Run Tests
```bash
make test
```

### Clean Temporary Files
```bash
make clean
```

## 📋 Command Line Options

| Option | Short | Description | Required |
|--------|-------|-------------|----------|
| `--genome` | `-g` | Path to genome FASTA file | Yes |
| `--annotation` | `-a` | Path to gene annotation GFF file | Yes |
| `--output` | `-o` | Output CSV file path | No |
| `--verbose` | `-v` | Enable verbose output | No |
| `--help` | `-h` | Show help message | No |

## 🧪 Test Data

The repository includes sample data for testing:
- `test_data/genome_sample.fasta`: Small E. coli genome fragment
- `test_data/annotation_sample.gff`: Corresponding gene annotations in GFF3 format

## 🐛 Error Handling

The tool validates input files and provides helpful error messages:

```bash
# If genome file doesn't exist
❌ Error: Genome file not found: missing_file.fasta

# If annotation file doesn't exist  
❌ Error: Annotation file not found: missing_annotation.gff

# If required arguments are missing
python genomnomnom.py
usage: genomnomnom.py [-h] --genome GENOME --annotation ANNOTATION [--output OUTPUT] [--verbose]
genomnomnom.py: error: the following arguments are required: --genome/-g, --annotation/-a
```

## 📊 Understanding the Output

### Genome Statistics
- **Total length**: Size of the genome in base pairs
- **GC content**: Percentage of G and C nucleotides
- **Number of contigs**: Number of contiguous sequences

### Start Codon Usage
- **ATG**: Most common start codon (methionine)
- **GTG**: Alternative start codon (valine)
- **TTG**: Less common start codon (leucine)

### Stop Codon Usage
- **TAA**: Amber stop codon
- **TAG**: Ochre stop codon  
- **TGA**: Opal stop codon

### ORF Statistics
- **Total ORFs**: Number of open reading frames found
- **Mean/Median length**: Average and middle ORF lengths
- **Longest/Shortest ORF**: Extreme ORF sizes

---

*Note: This is currently a mock implementation for demonstration purposes. Real genome analysis functionality will be added in future versions.*
