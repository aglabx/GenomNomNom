# GenomNomNom Usage Examples

## 🚀 Quick Start

### 1. Local File Analysis
```bash
python genomnomnom.py --genome test_data/genome_sample.fasta --annotation test_data/annotation_sample.gff
```

### 2. NCBI Genome Search and Download
```bash
python genomnomnom.py --species "Mycoplasma genitalium" --email your.email@example.com
```

Expected interactive output:
```
🧬 Starting GenomNomNom analysis...
========================================
🔍 NCBI mode: searching for 'Mycoplasma genitalium'
🔍 Searching NCBI for 'Mycoplasma genitalium' genomes...
✅ Found 1 reference genome(s)

📋 Available reference genomes:
--------------------------------------------------------------------------------
 1. Mycoplasmoides genitalium (bacteria)
    Accession: GCF_040556925.1
    Name: ASM4055692v1
    Level: N/A
    Size: 0.6 MB
    Date: 2024/10/16

Select genome (1-1) or 'q' to quit: 1

✅ Selected: Mycoplasmoides genitalium (bacteria) (GCF_040556925.1)
📥 Downloading FASTA file...
📥 Downloading GFF file...
✅ Files downloaded successfully!
📁 Using downloaded files:
   Genome: downloads/GCF_040556925.1_genomic.fna
   Annotation: downloads/GCF_040556925.1_genomic.gff

🔍 Parsing genome file: downloads/GCF_040556925.1_genomic.fna
✅ Genome parsed successfully
   📊 Found 1 sequence(s)
   📏 Total length: 580,076 bp
   🧬 GC content: 31.7%
```

### 3. Local Analysis - Expected Output
```
🧬 GenomNomNom Analysis Report
=====================================

📊 Genome Statistics:
- Total length: 623 bp
- GC content: 48.3%
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

### 4. Save Results to CSV
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

## 🌐 NCBI Integration Examples

### 1. Search for Bacterial Genomes
```bash
# Mycoplasma - one of the smallest bacterial genomes
python genomnomnom.py --species "Mycoplasma genitalium" --email your.email@example.com

# E. coli - model organism with well-annotated genome
python genomnomnom.py --species "Escherichia coli" --email your.email@example.com

# Bacillus subtilis - another model bacterium
python genomnomnom.py --species "Bacillus subtilis" --email your.email@example.com
```

### 2. Search for Eukaryotic Genomes
```bash
# Human reference genome (warning: very large!)
python genomnomnom.py --species "Homo sapiens" --email your.email@example.com

# Model organism yeast
python genomnomnom.py --species "Saccharomyces cerevisiae" --email your.email@example.com

# Fruit fly
python genomnomnom.py --species "Drosophila melanogaster" --email your.email@example.com
```

### 3. Save NCBI Results to CSV
```bash
python genomnomnom.py --species "Mycoplasma genitalium" --email your.email@example.com --output myco_analysis.csv
```

### 4. Test NCBI Functionality (Non-interactive)
```bash
python test_ncbi.py
```

Expected output:
```
🧪 Testing NCBI search functionality...
🔍 Searching for 'Mycoplasma genitalium'...
🔍 Searching NCBI for 'Mycoplasma genitalium' genomes...
✅ Found 1 reference genome(s)
✅ Found 1 assemblies

📋 Available assemblies:
 1. Mycoplasmoides genitalium (bacteria)
    Accession: GCF_040556925.1
    Name: ASM4055692v1
    Level: N/A
    Size: 0.6 MB

🔬 Testing download URLs construction...
   FASTA URL: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/556/925/GCF_040556925.1_ASM4055692v1/GCF_040556925.1_ASM4055692v1_genomic.fna.gz
   GFF URL: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/556/925/GCF_040556925.1_ASM4055692v1/GCF_040556925.1_ASM4055692v1_genomic.gff.gz
   FASTA file exists: True
   GFF file exists: True
✅ NCBI search test completed successfully!
```

### 5. NCBI Search Tips

**Good species for testing:**
- `"Mycoplasma genitalium"` - Very small genome (~580 kb)
- `"Escherichia coli"` - Standard model organism (~4.6 Mb)
- `"Bacillus subtilis"` - Gram-positive model (~4.2 Mb)

**Email requirement:**
- NCBI requires an email address for API access
- Use your real email address
- This helps NCBI contact you if there are issues with your queries

**Downloaded files location:**
- Files are saved in the `downloads/` directory
- FASTA files: `{accession}_genomic.fna`
- GFF files: `{accession}_genomic.gff`

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
