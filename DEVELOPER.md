# GenomNomNom Developer Documentation

## 📁 Project Structure

```
GenomNomNom/
├── genomnomnom.py          # Main application script
├── test_genomnomnom.py     # Test runner script
├── requirements.txt        # Python dependencies
├── Makefile               # Build automation
├── README.md              # User documentation
├── LICENSE                # License file
├── test_data/             # Sample data for testing
│   ├── genome_sample.fasta     # Sample genome FASTA
│   └── annotation_sample.gff   # Sample GFF annotation
└── drafts/
    └── paper.md           # Research paper draft
```

## 🏗️ Architecture

### Current Mock Implementation

The current version is a **mock implementation** that demonstrates the expected interface and output format. All analysis functions return hardcoded realistic data.

### Core Classes

1. **GenomeParser**: Handles FASTA genome file parsing
   - Currently validates file existence
   - Returns mock genome statistics

2. **AnnotationParser**: Handles GFF annotation file parsing  
   - Currently validates file existence
   - Returns mock gene counts

3. **CodonAnalyzer**: Analyzes codon usage patterns
   - Returns mock start/stop codon counts
   - Calculates realistic percentages

4. **ORFAnalyzer**: Analyzes Open Reading Frames
   - Returns mock ORF statistics

5. **ReportGenerator**: Creates output reports
   - Generates formatted console output
   - Exports results to CSV

## 🧪 Testing

### Running Tests
```bash
# Run all tests
python test_genomnomnom.py

# Or use Makefile
make test
```

### Test Data
- `test_data/genome_sample.fasta`: Small E. coli genome fragment
- `test_data/annotation_sample.gff`: Corresponding GFF3 annotations

## 🚀 Future Implementation Plan

### Phase 1: Real Parsing
- [ ] Implement actual FASTA parsing using BioPython
- [ ] Implement GFF3 parsing
- [ ] Calculate real genome statistics (length, GC content)

### Phase 2: Codon Analysis
- [ ] Extract coding sequences from genome + annotations
- [ ] Count actual start/stop codons
- [ ] Calculate codon usage statistics

### Phase 3: ORF Detection
- [ ] Implement ORF finder for all 6 reading frames
- [ ] Calculate ORF length distributions
- [ ] Compare annotated vs. predicted ORFs

### Phase 4: Advanced Features
- [ ] HMM-free gene prediction
- [ ] Comparative genomics across species
- [ ] Interactive visualizations
- [ ] Jupyter notebook integration

## 🛠️ Development Workflow

### Adding New Features

1. **Create branch**: `git checkout -b feature/new-feature`
2. **Implement**: Add functionality to appropriate class
3. **Update tests**: Add tests for new functionality
4. **Update docs**: Update README and this file
5. **Test**: Run `make test` to verify everything works
6. **Submit PR**: Create pull request for review

### Code Style

- Use descriptive variable names
- Add type hints where appropriate
- Include docstrings for all classes and methods
- Use emojis in user-facing output for fun! 🧬
- Follow PEP 8 style guidelines

### Mock Data Guidelines

Current mock data is based on realistic E. coli statistics:
- Genome length: ~4.6 Mbp
- GC content: ~50.8%  
- Start codon distribution: ATG (95.6%), GTG (3.8%), TTG (0.6%)
- Stop codon distribution: TAA (45.6%), TGA (47.2%), TAG (7.2%)

When updating mock data, ensure it remains biologically realistic.

## 🐛 Known Issues

1. Type annotations cause linter warnings (will be fixed in real implementation)
2. Pandas dependency not actually used yet (needed for future CSV handling)
3. Mock data is hardcoded (will be replaced with real analysis)

## 📝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Update documentation
6. Submit a pull request

For questions, please open an issue or contact the development team.

---

*Happy genome munching! 🍽️*
