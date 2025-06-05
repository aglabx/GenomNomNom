# GenomNomNom Makefile

.PHONY: help install test run clean

help:
	@echo "🧬 GenomNomNom - Available commands:"
	@echo "  make install    - Install dependencies"
	@echo "  make test       - Run test with sample data"
	@echo "  make run        - Run with sample data"
	@echo "  make clean      - Clean temporary files"
	@echo "  make help       - Show this help"

install:
	@echo "📦 Installing dependencies..."
	pip install -r requirements.txt

test:
	@echo "🧪 Running tests..."
	python test_genomnomnom.py

run:
	@echo "🚀 Running GenomNomNom with sample data..."
	python genomnomnom.py --genome test_data/genome_sample.fasta --annotation test_data/annotation_sample.gff --output results.csv --verbose

clean:
	@echo "🧹 Cleaning temporary files..."
	rm -f *.csv
	rm -rf __pycache__
	rm -rf .pytest_cache
	find . -name "*.pyc" -delete

demo:
	@echo "🎭 GenomNomNom Demo"
	@echo "=================="
	@echo "This will show the help first, then run analysis:"
	@echo ""
	python genomnomnom.py --help
	@echo ""
	@echo "Now running analysis with sample data..."
	@echo ""
	python genomnomnom.py --genome test_data/genome_sample.fasta --annotation test_data/annotation_sample.gff
