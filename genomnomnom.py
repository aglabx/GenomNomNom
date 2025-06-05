#!/usr/bin/env python3
"""
🧬 GenomNomNom - A playful tool for munching through genomes

This is the main entry point for the GenomNomNom toolkit.
Phase 1: Real FASTA and GFF parsing implementation.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import re
import csv

# Mock data for demonstration
MOCK_CODON_COUNTS = {
    'start_codons': {
        'ATG': 3847,
        'GTG': 152,
        'TTG': 24
    },
    'stop_codons': {
        'TAA': 1834,
        'TAG': 289,
        'TGA': 1900
    }
}

MOCK_ORF_STATS = {
    'total_orfs': 4023,
    'mean_length': 924.3,
    'median_length': 678,
    'longest_orf': 4932,
    'shortest_orf': 96
}

class GenomeParser:
    """Handles parsing of FASTA genome files"""
    
    def __init__(self, genome_file: str):
        self.genome_file = Path(genome_file)
        self.sequences = []
        self.stats = {}
    
    def parse(self) -> Dict:
        """Parse genome FASTA file using BioPython"""
        print(f"🔍 Parsing genome file: {self.genome_file}")
        
        # Check if file exists
        if not self.genome_file.exists():
            raise FileNotFoundError(f"Genome file not found: {self.genome_file}")
        
        try:
            # Parse FASTA file using BioPython
            self.sequences = list(SeqIO.parse(str(self.genome_file), "fasta"))
            
            if not self.sequences:
                raise ValueError(f"No sequences found in FASTA file: {self.genome_file}")
            
            # Calculate real statistics
            self.stats = self._calculate_stats()
            
            print("✅ Genome parsed successfully")
            print(f"   📊 Found {len(self.sequences)} sequence(s)")
            print(f"   📏 Total length: {self.stats['total_length']:,} bp")
            print(f"   🧬 GC content: {self.stats['gc_content']:.1f}%")
            
            return self.stats
            
        except Exception as e:
            raise ValueError(f"Error parsing FASTA file: {e}")
    
    def _calculate_stats(self) -> Dict:
        """Calculate genome statistics from parsed sequences"""
        total_length = 0
        total_gc_count = 0
        
        for seq_record in self.sequences:
            seq_str = str(seq_record.seq).upper()
            total_length += len(seq_str)
            
            # Count GC content
            gc_count = seq_str.count('G') + seq_str.count('C')
            total_gc_count += gc_count
        
        # Calculate GC percentage
        gc_content = (total_gc_count / total_length) * 100 if total_length > 0 else 0
        
        return {
            'total_length': total_length,
            'gc_content': round(gc_content, 1),
            'num_contigs': len(self.sequences),
            'sequences': self.sequences  # Store for later use
        }

class AnnotationParser:
    """Handles parsing of GFF annotation files"""
    
    def __init__(self, annotation_file: str):
        self.annotation_file = Path(annotation_file)
        self.genes = []
        self.cds_features = []
        self.stats = {}
    
    def parse(self) -> Dict:
        """Parse GFF annotation file"""
        print(f"🔍 Parsing annotation file: {self.annotation_file}")
        
        # Check if file exists
        if not self.annotation_file.exists():
            raise FileNotFoundError(f"Annotation file not found: {self.annotation_file}")
        
        try:
            # Simple GFF3 parser (basic implementation)
            self._parse_gff3()
            
            # Calculate statistics
            self.stats = self._calculate_stats()
            
            print("✅ Annotations parsed successfully")
            print(f"   🧬 Found {self.stats['total_genes']} genes")
            print(f"   💻 Found {self.stats['coding_genes']} coding sequences")
            
            return self.stats
            
        except Exception as e:
            raise ValueError(f"Error parsing GFF file: {e}")
    
    def _parse_gff3(self):
        """Parse GFF3 file format"""
        with open(self.annotation_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip comments and empty lines
                if line.startswith('#') or not line:
                    continue
                
                # Split GFF fields
                fields = line.split('\t')
                if len(fields) < 9:
                    continue
                
                seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
                
                # Store genes and CDS features
                if feature_type.lower() == 'gene':
                    self.genes.append({
                        'seqid': seqid,
                        'type': feature_type,
                        'start': int(start),
                        'end': int(end),
                        'strand': strand,
                        'attributes': attributes
                    })
                elif feature_type.lower() == 'cds':
                    self.cds_features.append({
                        'seqid': seqid,
                        'type': feature_type,
                        'start': int(start),
                        'end': int(end),
                        'strand': strand,
                        'phase': phase,
                        'attributes': attributes
                    })
    
    def _calculate_stats(self) -> Dict:
        """Calculate annotation statistics"""
        return {
            'total_genes': len(self.genes),
            'coding_genes': len(self.cds_features),
            'genes': self.genes,  # Store for later use
            'cds_features': self.cds_features  # Store for later use
        }

class CodonAnalyzer:
    """Analyzes codon usage patterns"""
    
    def __init__(self, genome_parser: GenomeParser, annotation_parser: AnnotationParser):
        self.genome_parser = genome_parser
        self.annotation_parser = annotation_parser
        self.codon_counts = {}
    
    def count_codons(self) -> Dict:
        """Count start and stop codons (MOCK implementation)"""
        print("🔬 Analyzing codon usage...")
        
        # Mock implementation
        self.codon_counts = MOCK_CODON_COUNTS
        
        print("✅ Codon analysis complete")
        return self.codon_counts
    
    def calculate_percentages(self) -> Dict:
        """Calculate codon usage percentages"""
        if not self.codon_counts:
            self.count_codons()
        
        results = {}
        
        # Calculate start codon percentages
        total_start = sum(self.codon_counts['start_codons'].values())
        results['start_percentages'] = {
            codon: (count / total_start) * 100 
            for codon, count in self.codon_counts['start_codons'].items()
        }
        
        # Calculate stop codon percentages
        total_stop = sum(self.codon_counts['stop_codons'].values())
        results['stop_percentages'] = {
            codon: (count / total_stop) * 100 
            for codon, count in self.codon_counts['stop_codons'].items()
        }
        
        return results

class ORFAnalyzer:
    """Analyzes Open Reading Frames"""
    
    def __init__(self, genome_parser: GenomeParser, annotation_parser: AnnotationParser):
        self.genome_parser = genome_parser
        self.annotation_parser = annotation_parser
        self.orf_stats = {}
    
    def analyze_orfs(self) -> Dict:
        """Analyze ORF statistics (MOCK implementation)"""
        print("🔬 Analyzing ORFs...")
        
        # Mock implementation
        self.orf_stats = MOCK_ORF_STATS
        
        print("✅ ORF analysis complete")
        return self.orf_stats

class ReportGenerator:
    """Generates output reports"""
    
    def __init__(self):
        self.results = {}
    
    def generate_summary(self, genome_stats: Dict, codon_stats: Dict, orf_stats: Dict) -> str:
        """Generate summary report"""
        
        codon_analyzer = CodonAnalyzer(None, None)
        codon_analyzer.codon_counts = codon_stats
        percentages = codon_analyzer.calculate_percentages()
        
        summary = f"""
🧬 GenomNomNom Analysis Report
=====================================

📊 Genome Statistics:
- Total length: {genome_stats['total_length']:,} bp
- GC content: {genome_stats['gc_content']:.1f}%
- Number of contigs: {genome_stats['num_contigs']}

🚀 Start Codon Usage:
- ATG: {codon_stats['start_codons']['ATG']:,} ({percentages['start_percentages']['ATG']:.1f}%)
- GTG: {codon_stats['start_codons']['GTG']:,} ({percentages['start_percentages']['GTG']:.1f}%)
- TTG: {codon_stats['start_codons']['TTG']:,} ({percentages['start_percentages']['TTG']:.1f}%)

🛑 Stop Codon Usage:
- TAA: {codon_stats['stop_codons']['TAA']:,} ({percentages['stop_percentages']['TAA']:.1f}%)
- TAG: {codon_stats['stop_codons']['TAG']:,} ({percentages['stop_percentages']['TAG']:.1f}%)
- TGA: {codon_stats['stop_codons']['TGA']:,} ({percentages['stop_percentages']['TGA']:.1f}%)

📏 ORF Statistics:
- Total ORFs: {orf_stats['total_orfs']:,}
- Mean length: {orf_stats['mean_length']:.1f} bp
- Median length: {orf_stats['median_length']} bp
- Longest ORF: {orf_stats['longest_orf']:,} bp
- Shortest ORF: {orf_stats['shortest_orf']} bp

🍽️ Nom nom nom! Analysis complete.
"""
        return summary
    
    def save_csv(self, filename: str, codon_stats: Dict, orf_stats: Dict):
        """Save results to CSV file using csv module"""
        with open(filename, 'w', newline='') as csvfile:
            fieldnames = ['metric_type', 'codon', 'count']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            # Write header
            writer.writeheader()
            
            # Add codon counts
            for codon_type in ['start_codons', 'stop_codons']:
                for codon, count in codon_stats[codon_type].items():
                    writer.writerow({
                        'metric_type': codon_type,
                        'codon': codon,
                        'count': count
                    })
            
            # Add ORF stats
            for metric, value in orf_stats.items():
                writer.writerow({
                    'metric_type': 'orf_stats',
                    'codon': metric,
                    'count': value
                })
        
        print(f"📄 Results saved to: {filename}")

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="🧬 GenomNomNom - A playful tool for munching through genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python genomnomnom.py --genome genome.fasta --annotation genes.gff
  python genomnomnom.py --genome genome.fasta --annotation genes.gff --output results.csv
        """
    )
    
    parser.add_argument(
        "--genome", "-g",
        required=True,
        help="Path to genome FASTA file"
    )
    
    parser.add_argument(
        "--annotation", "-a", 
        required=True,
        help="Path to gene annotation file (GFF format)"
    )
    
    parser.add_argument(
        "--output", "-o",
        help="Output CSV file path (optional)"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose output"
    )
    
    args = parser.parse_args()
    
    try:
        print("🧬 Starting GenomNomNom analysis...")
        print("=" * 40)
        
        # Parse genome
        genome_parser = GenomeParser(args.genome)
        genome_stats = genome_parser.parse()
        
        # Parse annotations
        annotation_parser = AnnotationParser(args.annotation)
        annotation_stats = annotation_parser.parse()
        
        # Analyze codons
        codon_analyzer = CodonAnalyzer(genome_parser, annotation_parser)
        codon_stats = codon_analyzer.count_codons()
        
        # Analyze ORFs
        orf_analyzer = ORFAnalyzer(genome_parser, annotation_parser)
        orf_stats = orf_analyzer.analyze_orfs()
        
        # Generate report
        report_generator = ReportGenerator()
        summary = report_generator.generate_summary(genome_stats, codon_stats, orf_stats)
        
        print(summary)
        
        # Save CSV if requested
        if args.output:
            report_generator.save_csv(args.output, codon_stats, orf_stats)
        
        print("🎉 GenomNomNom finished successfully!")
        
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
