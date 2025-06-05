#!/usr/bin/env python3
"""
🧬 GenomNomNom - A playful tool for munching through genomes

This is the main entry point for the GenomNomNom toolkit.
Phase 1: Real FASTA and GFF parsing implementation.
Phase 2: NCBI genome search and download functionality.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import re
import csv
import requests
import json
import time
import os
from urllib.parse import quote
import gzip
import shutil

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

class NCBISearcher:
    """Handles searching and downloading genomes from NCBI"""
    
    def __init__(self, email: str = "your.email@example.com"):
        self.email = email
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.ftp_base = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
        
    def search_genomes(self, species_name: str, max_results: int = 10) -> List[Dict]:
        """Search for reference genomes by species name"""
        print(f"🔍 Searching NCBI for '{species_name}' genomes...")
        
        # Search for assembly records
        search_term = f'"{species_name}"[Organism] AND "reference genome"[RefSeq Category]'
        
        try:
            # Step 1: Search assembly database
            search_url = f"{self.base_url}/esearch.fcgi"
            search_params = {
                'db': 'assembly',
                'term': search_term,
                'retmax': max_results,
                'retmode': 'json',
                'email': self.email
            }
            
            response = requests.get(search_url, params=search_params)
            response.raise_for_status()
            search_results = response.json()
            
            if not search_results.get('esearchresult', {}).get('idlist'):
                print(f"❌ No reference genomes found for '{species_name}'")
                return []
            
            # Step 2: Get assembly details
            assembly_ids = search_results['esearchresult']['idlist']
            summary_url = f"{self.base_url}/esummary.fcgi"
            summary_params = {
                'db': 'assembly',
                'id': ','.join(assembly_ids),
                'retmode': 'json',
                'email': self.email
            }
            
            response = requests.get(summary_url, params=summary_params)
            response.raise_for_status()
            summary_results = response.json()
            
            # Parse assembly information
            assemblies = []
            for assembly_id in assembly_ids:
                if assembly_id in summary_results['result']:
                    assembly_info = summary_results['result'][assembly_id]
                    assemblies.append({
                        'id': assembly_id,
                        'accession': assembly_info.get('assemblyaccession', 'N/A'),
                        'name': assembly_info.get('assemblyname', 'N/A'),
                        'organism': assembly_info.get('organism', 'N/A'),
                        'level': assembly_info.get('assemblylevel', 'N/A'),
                        'size': assembly_info.get('totalsize', 0),
                        'submission_date': assembly_info.get('submissiondate', 'N/A'),
                        'ftp_path': assembly_info.get('ftppath_refseq', assembly_info.get('ftppath_genbank', ''))
                    })
            
            print(f"✅ Found {len(assemblies)} reference genome(s)")
            return assemblies
            
        except requests.RequestException as e:
            print(f"❌ Error searching NCBI: {e}")
            return []
        except Exception as e:
            print(f"❌ Unexpected error: {e}")
            return []
    
    def download_genome_files(self, assembly_info: Dict, output_dir: str = "downloads") -> Tuple[Optional[str], Optional[str]]:
        """Download FASTA and GFF files for a specific assembly"""
        if not assembly_info.get('ftp_path'):
            print("❌ No FTP path available for this assembly")
            return None, None
        
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        accession = assembly_info['accession']
        ftp_base = assembly_info['ftp_path'].replace('ftp://', 'https://')
        
        # Construct file URLs
        base_name = ftp_base.split('/')[-1]
        fasta_url = f"{ftp_base}/{base_name}_genomic.fna.gz"
        gff_url = f"{ftp_base}/{base_name}_genomic.gff.gz"
        
        fasta_file = None
        gff_file = None
        
        try:
            # Download FASTA file
            print(f"📥 Downloading FASTA file...")
            fasta_file = self._download_and_extract(fasta_url, output_path, f"{accession}_genomic.fna")
            
            # Download GFF file
            print(f"📥 Downloading GFF file...")
            gff_file = self._download_and_extract(gff_url, output_path, f"{accession}_genomic.gff")
            
            print("✅ Files downloaded successfully!")
            return fasta_file, gff_file
            
        except Exception as e:
            print(f"❌ Error downloading files: {e}")
            return None, None
    
    def _download_and_extract(self, url: str, output_dir: Path, filename: str) -> str:
        """Download and extract a gzipped file"""
        compressed_file = output_dir / f"{filename}.gz"
        extracted_file = output_dir / filename
        
        # Download compressed file
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        with open(compressed_file, 'wb') as f:
            shutil.copyfileobj(response.raw, f)
        
        # Extract file
        with gzip.open(compressed_file, 'rb') as f_in:
            with open(extracted_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        # Remove compressed file
        compressed_file.unlink()
        
        return str(extracted_file)
    
    def interactive_genome_selection(self, species_name: str) -> Tuple[Optional[str], Optional[str]]:
        """Interactive workflow for selecting and downloading a genome"""
        assemblies = self.search_genomes(species_name)
        
        if not assemblies:
            return None, None
        
        print("\n📋 Available reference genomes:")
        print("-" * 80)
        for i, assembly in enumerate(assemblies, 1):
            size_mb = int(assembly['size']) / 1_000_000 if assembly['size'] else 0
            print(f"{i:2d}. {assembly['organism']}")
            print(f"    Accession: {assembly['accession']}")
            print(f"    Name: {assembly['name']}")
            print(f"    Level: {assembly['level']}")
            print(f"    Size: {size_mb:.1f} MB")
            print(f"    Date: {assembly['submission_date']}")
            print()
        
        # User selection
        while True:
            try:
                choice = input(f"Select genome (1-{len(assemblies)}) or 'q' to quit: ").strip()
                if choice.lower() == 'q':
                    return None, None
                
                choice_idx = int(choice) - 1
                if 0 <= choice_idx < len(assemblies):
                    selected_assembly = assemblies[choice_idx]
                    print(f"\n✅ Selected: {selected_assembly['organism']} ({selected_assembly['accession']})")
                    
                    # Download files
                    fasta_file, gff_file = self.download_genome_files(selected_assembly)
                    return fasta_file, gff_file
                else:
                    print("❌ Invalid selection. Please try again.")
            except ValueError:
                print("❌ Please enter a number or 'q' to quit.")
            except KeyboardInterrupt:
                print("\n👋 Cancelled by user")
                return None, None

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
  # Analyze local files
  python genomnomnom.py --genome genome.fasta --annotation genes.gff
  python genomnomnom.py --genome genome.fasta --annotation genes.gff --output results.csv
  
  # Search and download from NCBI
  python genomnomnom.py --species "Escherichia coli" --email your.email@example.com
  python genomnomnom.py --species "Homo sapiens" --email your.email@example.com --output human_results.csv
        """
    )
    
    # Local file analysis arguments
    parser.add_argument(
        "--genome", "-g",
        help="Path to genome FASTA file"
    )
    
    parser.add_argument(
        "--annotation", "-a", 
        help="Path to gene annotation file (GFF format) - required for local mode"
    )
    
    # NCBI search arguments
    parser.add_argument(
        "--species", "-s",
        help="Species name to search for (e.g., 'Escherichia coli')"
    )
    
    parser.add_argument(
        "--email", "-e",
        help="Email address for NCBI API (required for species search)"
    )
    
    # Common arguments
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
    
    # Validate arguments based on mode
    if args.genome and not args.annotation:
        parser.error("--annotation is required when using --genome")
    
    if args.species and not args.email:
        parser.error("--email is required when using --species")
    
    if not args.genome and not args.species:
        parser.error("Either --genome or --species must be specified")
    
    if args.genome and args.species:
        parser.error("Cannot specify both --genome and --species. Choose one mode.")
    
    try:
        print("🧬 Starting GenomNomNom analysis...")
        print("=" * 40)
        
        genome_file = args.genome
        annotation_file = args.annotation
        
        # NCBI search and download mode
        if args.species:
            print(f"🔍 NCBI mode: searching for '{args.species}'")
            
            ncbi_searcher = NCBISearcher(email=args.email)
            genome_file, annotation_file = ncbi_searcher.interactive_genome_selection(args.species)
            
            if not genome_file or not annotation_file:
                print("❌ No genome selected or download failed. Exiting.")
                sys.exit(1)
            
            print(f"📁 Using downloaded files:")
            print(f"   Genome: {genome_file}")
            print(f"   Annotation: {annotation_file}")
            print()
        
        # Parse genome
        genome_parser = GenomeParser(genome_file)
        genome_stats = genome_parser.parse()
        
        # Parse annotations  
        annotation_parser = AnnotationParser(annotation_file)
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
