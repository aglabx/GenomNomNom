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

# Note: Mock data removed - now using real analysis only

# Genetic code table (Standard genetic code)
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# Amino acid names
AMINO_ACID_NAMES = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': 'Stop'
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
            
            # Debug: Print available fields for the first assembly
            if assembly_ids and assembly_ids[0] in summary_results['result']:
                first_assembly = summary_results['result'][assembly_ids[0]]
                print("🔧 Debug - Available assembly fields:")
                for key in sorted(first_assembly.keys()):
                    value = first_assembly[key]
                    print(f"   {key}: {value}")
                print()
            
            # Parse assembly information
            assemblies = []
            for assembly_id in assembly_ids:
                if assembly_id in summary_results['result']:
                    assembly_info = summary_results['result'][assembly_id]
                    
                    # Try different field names for assembly level
                    level = (assembly_info.get('assemblylevel') or 
                            assembly_info.get('assembly_level') or 
                            assembly_info.get('level') or
                            assembly_info.get('assemblytype') or
                            assembly_info.get('assemblyclass') or
                            assembly_info.get('assemblyDescription', 'N/A'))
                    
                    # Try different field names for genome size
                    size = (assembly_info.get('totalsize') or 
                           assembly_info.get('total_length') or 
                           assembly_info.get('genome_size') or 
                           assembly_info.get('totallength') or
                           assembly_info.get('assemblysize') or
                           assembly_info.get('stat_total_sequence_length') or
                           0)
                    
                    assemblies.append({
                        'id': assembly_id,
                        'accession': assembly_info.get('assemblyaccession', 'N/A'),
                        'name': assembly_info.get('assemblyname', 'N/A'),
                        'organism': assembly_info.get('organism', 'N/A'),
                        'level': level,
                        'size': size,
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
        self.all_codon_counts = {}
        self.start_codons = ['ATG', 'GTG', 'TTG']
        self.stop_codons = ['TAA', 'TAG', 'TGA']
    
    def count_codons(self) -> Dict:
        """Count start and stop codons based on real CDS sequences"""
        print("🔬 Analyzing codon usage...")
        
        start_counts = {codon: 0 for codon in self.start_codons}
        stop_counts = {codon: 0 for codon in self.stop_codons}
        
        # Get CDS features from annotation
        cds_features = self.annotation_parser.stats.get('cds_features', [])
        
        if not cds_features:
            print("⚠️  No CDS features found in annotation file")
            print("❌ Cannot perform codon analysis without CDS annotations")
            # Return empty counts instead of mock data
            self.codon_counts = {
                'start_codons': {'ATG': 0, 'GTG': 0, 'TTG': 0},
                'stop_codons': {'TAA': 0, 'TAG': 0, 'TGA': 0}
            }
            return self.codon_counts
        
        # Get genome sequences (as dictionary for quick lookup)
        sequences = {}
        for seq_record in self.genome_parser.sequences:
            sequences[seq_record.id] = str(seq_record.seq).upper()
        
        cds_analyzed = 0
        for cds in cds_features:
            seqid = cds['seqid']
            start = cds['start'] - 1  # Convert to 0-based indexing
            end = cds['end']
            strand = cds['strand']
            
            if seqid not in sequences:
                continue
            
            # Extract CDS sequence
            cds_seq = sequences[seqid][start:end]
            
            # Handle reverse strand
            if strand == '-':
                cds_seq = self._reverse_complement(cds_seq)
            
            # Skip if sequence is too short
            if len(cds_seq) < 6:  # Need at least start + stop codon
                continue
            
            # Count start codon (first 3 bases)
            start_codon = cds_seq[:3]
            if start_codon in start_counts:
                start_counts[start_codon] += 1
            
            # Count stop codon (last 3 bases)
            if len(cds_seq) >= 3:
                stop_codon = cds_seq[-3:]
                if stop_codon in stop_counts:
                    stop_counts[stop_codon] += 1
            
            cds_analyzed += 1
        
        self.codon_counts = {
            'start_codons': start_counts,
            'stop_codons': stop_counts
        }
        
        print(f"   📊 Analyzed {cds_analyzed} CDS sequences")
        print("✅ Codon analysis complete")
        return self.codon_counts
    
    def _reverse_complement(self, seq: str) -> str:
        """Calculate reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(seq))
    
    def calculate_percentages(self) -> Dict:
        """Calculate codon usage percentages"""
        if not self.codon_counts:
            self.count_codons()
        
        results = {}
        
        # Calculate start codon percentages
        total_start = sum(self.codon_counts['start_codons'].values())
        if total_start > 0:
            results['start_percentages'] = {
                codon: (count / total_start) * 100 
                for codon, count in self.codon_counts['start_codons'].items()
            }
        else:
            results['start_percentages'] = {codon: 0 for codon in self.start_codons}
        
        # Calculate stop codon percentages
        total_stop = sum(self.codon_counts['stop_codons'].values())
        if total_stop > 0:
            results['stop_percentages'] = {
                codon: (count / total_stop) * 100 
                for codon, count in self.codon_counts['stop_codons'].items()
            }
        else:
            results['stop_percentages'] = {codon: 0 for codon in self.stop_codons}
        
        return results

    def analyze_all_codons(self) -> None:
        """Analyze all 64 codons from CDS sequences"""
        print("🔬 Analyzing all 64 codons...")
        
        # Initialize counts for all possible codons
        all_codons = []
        bases = ['A', 'T', 'G', 'C']
        for base1 in bases:
            for base2 in bases:
                for base3 in bases:
                    all_codons.append(base1 + base2 + base3)
        
        self.all_codon_counts = {codon: 0 for codon in all_codons}
        
        # Get CDS features from annotation
        cds_features = self.annotation_parser.stats.get('cds_features', [])
        
        if not cds_features:
            print("⚠️  No CDS features found for all-codon analysis")
            return
        
        # Get genome sequences
        sequences = {}
        for seq_record in self.genome_parser.sequences:
            sequences[seq_record.id] = str(seq_record.seq).upper()
        
        total_codons_found = 0
        cds_analyzed = 0
        
        for cds in cds_features:
            seqid = cds['seqid']
            start = cds['start'] - 1  # Convert to 0-based indexing
            end = cds['end']
            strand = cds['strand']
            
            if seqid not in sequences:
                continue
            
            # Extract CDS sequence
            cds_seq = sequences[seqid][start:end]
            
            # Handle reverse strand
            if strand == '-':
                cds_seq = self._reverse_complement(cds_seq)
            
            # Count all codons in this CDS (in triplets)
            for i in range(0, len(cds_seq) - 2, 3):
                codon = cds_seq[i:i+3]
                if len(codon) == 3 and all(base in 'ATGC' for base in codon):
                    if codon in self.all_codon_counts:
                        self.all_codon_counts[codon] += 1
                        total_codons_found += 1
            
            cds_analyzed += 1
        
        print(f"   📊 Analyzed {cds_analyzed} CDS sequences")
        print(f"   📊 Found {total_codons_found:,} total codons")
        print("✅ All-codon analysis complete")

    def generate_codon_table_report(self) -> str:
        """Generate a beautiful codon usage table report"""
        if not self.all_codon_counts:
            self.analyze_all_codons()
        
        total_codons = sum(self.all_codon_counts.values())
        if total_codons == 0:
            return "⚠️ No codon data available for table generation"
        
        # Group codons by amino acid
        aa_groups = {}
        for codon, aa in GENETIC_CODE.items():
            if aa not in aa_groups:
                aa_groups[aa] = []
            aa_groups[aa].append(codon)
        
        # Calculate frequencies and identify unusual patterns
        codon_frequencies = {}
        aa_frequencies = {}
        unusual_codons = []
        rare_codons = []
        dominant_codons = []
        
        for codon, count in self.all_codon_counts.items():
            freq = (count / total_codons) * 100 if total_codons > 0 else 0
            codon_frequencies[codon] = freq
            
            aa = GENETIC_CODE[codon]
            if aa not in aa_frequencies:
                aa_frequencies[aa] = 0
            aa_frequencies[aa] += freq
        
        # Identify unusual codon usage patterns
        for aa, codons in aa_groups.items():
            if len(codons) > 1:  # Amino acids with multiple codons
                aa_total = sum(codon_frequencies[c] for c in codons)
                if aa_total > 0:
                    for codon in codons:
                        relative_usage = (codon_frequencies[codon] / aa_total) * 100
                        # Flag as unusual patterns
                        if relative_usage < 3 and codon_frequencies[codon] > 0.05:
                            rare_codons.append((codon, aa, relative_usage, codon_frequencies[codon]))
                        elif relative_usage > 85:
                            dominant_codons.append((codon, aa, relative_usage, codon_frequencies[codon]))
                        elif (relative_usage < 5 and codon_frequencies[codon] > 0.1) or relative_usage > 80:
                            unusual_codons.append((codon, aa, relative_usage, codon_frequencies[codon]))
        
        # Build the report
        report = "\n" + "🧬" * 20 + " CODON USAGE ANALYSIS " + "🧬" * 20 + "\n"
        report += "=" * 80 + "\n"
        
        # Sort amino acids by frequency (excluding stop codons for main table)
        sorted_aas = sorted([(aa, freq) for aa, freq in aa_frequencies.items() if aa != '*'], 
                           key=lambda x: x[1], reverse=True)
        
        report += "\n📊 CODON USAGE BY AMINO ACID (sorted by frequency)\n"
        report += "─" * 80 + "\n"
        
        for i, (aa, aa_freq) in enumerate(sorted_aas):
            aa_name = AMINO_ACID_NAMES.get(aa, aa)
            
            # Add ranking emojis
            rank_emoji = ""
            if i == 0:
                rank_emoji = "🥇 "
            elif i == 1:
                rank_emoji = "🥈 "
            elif i == 2:
                rank_emoji = "🥉 "
            
            # Color coding based on frequency
            freq_emoji = ""
            if aa_freq > 7:
                freq_emoji = "🔴"  # Very high
            elif aa_freq > 5:
                freq_emoji = "🟠"  # High
            elif aa_freq > 3:
                freq_emoji = "🟡"  # Medium
            elif aa_freq > 1:
                freq_emoji = "🟢"  # Low
            else:
                freq_emoji = "🔵"  # Very low
                
            report += f"\n{rank_emoji}{freq_emoji} {aa_name} ({aa}) — {aa_freq:.2f}% of all codons\n"
            
            # Get codons for this amino acid, sorted by frequency
            aa_codons = [(c, codon_frequencies[c]) for c in aa_groups[aa]]
            aa_codons.sort(key=lambda x: x[1], reverse=True)
            
            # Create a nice table for codons
            report += "   ┌" + "─" * 70 + "┐\n"
            
            for j, (codon, freq) in enumerate(aa_codons):
                count = self.all_codon_counts[codon]
                # Calculate relative usage within this amino acid
                relative = (freq / aa_freq) * 100 if aa_freq > 0 else 0
                
                # Add visual indicators
                indicators = ""
                if freq > 4.0:  # Very high usage
                    indicators += "🔥"
                if freq < 0.05:  # Very low usage
                    indicators += "💤"
                if relative > 85:  # Extremely dominant codon for this AA
                    indicators += "�"
                elif relative > 70:  # Dominant codon for this AA
                    indicators += "⭐"
                if relative < 3 and freq > 0.05:  # Rare but present
                    indicators += "🔸"
                
                # Create progress bar for relative usage
                bar_length = 20
                filled = int((relative / 100) * bar_length)
                bar = "█" * filled + "░" * (bar_length - filled)
                
                report += f"   │ {codon}: {count:>6,} ({freq:>5.2f}% total, {relative:>5.1f}% of {aa}) [{bar}] {indicators}\n"
            
            report += "   └" + "─" * 70 + "┘\n"
        
        # Stop codons section with special formatting
        report += f"\n🛑 STOP CODONS — End of Translation\n"
        report += "─" * 40 + "\n"
        stop_codons = [(c, codon_frequencies[c]) for c in aa_groups['*']]
        stop_codons.sort(key=lambda x: x[1], reverse=True)
        
        report += "┌" + "─" * 38 + "┐\n"
        for codon, freq in stop_codons:
            count = self.all_codon_counts[codon]
            stop_emoji = "🔴" if freq > 1.5 else "🟡" if freq > 0.5 else "🟢"
            report += f"│ {stop_emoji} {codon}: {count:>6,} ({freq:>5.2f}% of all codons) │\n"
        report += "└" + "─" * 38 + "┘\n"
        
        # Enhanced unusual patterns section
        all_unusual = rare_codons + dominant_codons + unusual_codons
        if all_unusual:
            report += "\n⚠️  UNUSUAL CODON USAGE PATTERNS\n"
            report += "─" * 50 + "\n"
            
            if dominant_codons:
                report += "👑 DOMINANT CODONS (>85% usage for their amino acid):\n"
                for codon, aa, relative_usage, total_freq in sorted(dominant_codons, key=lambda x: x[2], reverse=True):
                    aa_name = AMINO_ACID_NAMES.get(aa, aa)
                    report += f"   {codon} ({aa_name}): {relative_usage:.1f}% dominance, {total_freq:.2f}% total\n"
                report += "\n"
            
            if rare_codons:
                report += "🔸 RARE CODONS (<3% usage for their amino acid):\n"
                for codon, aa, relative_usage, total_freq in sorted(rare_codons, key=lambda x: x[3], reverse=True):
                    aa_name = AMINO_ACID_NAMES.get(aa, aa)
                    report += f"   {codon} ({aa_name}): {relative_usage:.1f}% of {aa} codons, {total_freq:.2f}% total\n"
                report += "\n"
        
        # Enhanced summary with additional statistics
        report += "📈 COMPREHENSIVE SUMMARY\n"
        report += "─" * 30 + "\n"
        report += f"🔢 Total codons analyzed: {total_codons:,}\n"
        report += f"🥇 Most frequent AA: {sorted_aas[0][0]} ({AMINO_ACID_NAMES[sorted_aas[0][0]]}) — {sorted_aas[0][1]:.2f}%\n"
        report += f"🏁 Least frequent AA: {sorted_aas[-1][0]} ({AMINO_ACID_NAMES[sorted_aas[-1][0]]}) — {sorted_aas[-1][1]:.2f}%\n"
        
        # Calculate codon bias metrics
        max_codon_freq = max(codon_frequencies.values())
        min_codon_freq = min([f for f in codon_frequencies.values() if f > 0])
        
        report += f"📊 Codon frequency range: {min_codon_freq:.3f}% — {max_codon_freq:.2f}%\n"
        report += f"🎯 Codon usage bias ratio: {max_codon_freq/min_codon_freq:.1f}:1\n"
        
        # Count how many codons are actually used
        used_codons = sum(1 for freq in codon_frequencies.values() if freq > 0)
        report += f"✅ Codons in use: {used_codons}/64 ({used_codons/64*100:.1f}%)\n"
        
        return report
    
    def generate_classic_codon_table(self) -> str:
        """Generate a classic textbook-style square codon usage table with color visualization"""
        if not self.all_codon_counts:
            self.analyze_all_codons()
        
        total_codons = sum(self.all_codon_counts.values())
        if total_codons == 0:
            return "⚠️ No codon data available for classic table generation"
        
        # Calculate percentages for each codon
        codon_percentages = {}
        for codon, count in self.all_codon_counts.items():
            codon_percentages[codon] = (count / total_codons) * 100 if total_codons > 0 else 0
        
        # Color function for percentages with more contrast
        def colorize_percentage(percentage):
            """Apply color based on percentage value with high contrast"""
            if percentage >= 2.5:
                return f"\033[1;92m{percentage:4.1f}\033[0m"  # Bold bright green for high usage
            elif percentage >= 1.5:
                return f"\033[1;93m{percentage:4.1f}\033[0m"  # Bold bright yellow for medium usage
            elif percentage >= 0.5:
                return f"\033[1;96m{percentage:4.1f}\033[0m"  # Bold bright cyan for low usage
            elif percentage > 0:
                return f"\033[1;97m{percentage:4.1f}\033[0m"  # Bold bright white for very low usage
            else:
                return f"\033[1;90m{percentage:4.1f}\033[0m"  # Bold dark gray for zero usage
        
        bases = ['T', 'C', 'A', 'G']  # Using T instead of U for DNA
        
        report = "\n" + "🧬" * 30 + "\n"
        report += "📚 CLASSIC CODON USAGE TABLE (Textbook Style)\n"
        report += "🧬" * 30 + "\n"
        report += f"🔢 Total codons analyzed: \033[96m{total_codons:,}\033[0m\n"
        report += "═" * 78 + "\n\n"
        
        # Enhanced header with proper spacing to match data rows
        report += "       ┃"
        for second_base in bases:
            report += f"        \033[1;94m{second_base}\033[0m        ┃"
        report += "\n"
        
        report += "  \033[1;94m1st\033[0m  ┃"
        for _ in bases:
            report += "  \033[1;94m2nd position\033[0m        ┃"
        report += "\n"
        
        report += " \033[1;94mpos.\033[0m  ┃"
        for _ in bases:
            report += "  \033[1;94mT     C     A     G\033[0m  ┃"
        report += "\n"
        
        # Thicker separator matching the data width
        report += "━━━━━━━╋"
        for _ in bases:
            report += "━━━━━━━━━━━━━━━━━━━━━━━╋"
        report += "\n"
        
        # Data rows with enhanced formatting
        for first_base in bases:
            # Row with colored percentages
            report += f"   \033[1;94m{first_base}\033[0m   ┃"
            
            for second_base in bases:
                cell_data = ""
                for third_base in bases:
                    codon = first_base + second_base + third_base
                    percentage = codon_percentages.get(codon, 0)
                    cell_data += f"  {colorize_percentage(percentage)}"
                report += cell_data + " ┃"
            report += "\n"
            
            # Row with amino acids
            report += "       ┃"
            for second_base in bases:
                cell_data = ""
                for third_base in bases:
                    codon = first_base + second_base + third_base
                    aa = GENETIC_CODE.get(codon, '?')
                    if aa == '*':
                        aa_colored = f"\033[1;91mSTP\033[0m"  # Bold red for stop codons
                    else:
                        aa_colored = f"\033[1;37m{aa:>3}\033[0m"  # Bold white for amino acids
                    cell_data += f"  {aa_colored:5}"
                report += cell_data + " ┃"
            report += "\n"
            
            # Separator line between rows
            if first_base != bases[-1]:
                report += "       ┃"
                for _ in bases:
                    report += "                       ┃"
                report += "\n"
            else:
                # Final separator
                report += "━━━━━━━╋"
                for _ in bases:
                    report += "━━━━━━━━━━━━━━━━━━━━━━━╋"
                report += "\n"
        
        # Enhanced legend with more contrasting colors
        report += "\n🎨 \033[1;94mCOLOR LEGEND:\033[0m\n"
        report += "━" * 40 + "\n"
        report += f"• \033[1;92m■\033[0m High usage (≥2.5%)\n"
        report += f"• \033[1;93m■\033[0m Medium usage (1.5-2.4%)\n"
        report += f"• \033[1;96m■\033[0m Low usage (0.5-1.4%)\n"
        report += f"• \033[1;97m■\033[0m Very low usage (<0.5%)\n"
        report += f"• \033[1;90m■\033[0m Unused (0%)\n"
        report += f"• \033[1;91mSTP\033[0m Stop codons\n\n"
        
        report += "📋 \033[1;94mTABLE ORGANIZATION:\033[0m\n"
        report += "━" * 40 + "\n"
        report += "• Rows: 1st nucleotide position (T, C, A, G)\n"
        report += "• Columns: 2nd nucleotide position (T, C, A, G)\n"
        report += "• Within cells: 3rd position from left to right (T, C, A, G)\n"
        report += "• Numbers show percentage of total codon usage\n"
        report += "• Letters show corresponding amino acid codes\n"
        
        # Enhanced statistics
        sorted_codons = sorted(codon_percentages.items(), key=lambda x: x[1], reverse=True)
        most_used = sorted_codons[0] if sorted_codons else ('---', 0)
        used_codons = [c for c in sorted_codons if c[1] > 0]
        least_used = used_codons[-1] if used_codons else ('---', 0)
        unused_count = sum(1 for _, pct in sorted_codons if pct == 0)
        
        report += f"\n📊 \033[94mUSAGE STATISTICS:\033[0m\n"
        report += "━" * 40 + "\n"
        report += f"• Most used codon: \033[92m{most_used[0]}\033[0m (\033[92m{most_used[1]:.2f}%\033[0m)\n"
        if used_codons:
            report += f"• Least used codon: \033[37m{least_used[0]}\033[0m (\033[37m{least_used[1]:.2f}%\033[0m)\n"
        report += f"• Used codons: \033[96m{len(used_codons)}/64\033[0m\n"
        report += f"• Unused codons: \033[90m{unused_count}/64\033[0m\n"
        
        # Calculate codon usage evenness
        if used_codons:
            used_percentages = [pct for _, pct in used_codons]
            evenness = min(used_percentages) / max(used_percentages) if max(used_percentages) > 0 else 0
            if evenness >= 0.5:
                evenness_color = "\033[92m"  # Green
            elif evenness >= 0.2:
                evenness_color = "\033[93m"  # Yellow
            else:
                evenness_color = "\033[91m"  # Red
            report += f"• Usage evenness: {evenness_color}{evenness:.3f}\033[0m (1.0 = perfectly even)\n"
        
        return report

    def analyze_all_codons(self) -> None:
        """Analyze all 64 codons from CDS sequences"""
        print("🔬 Analyzing all 64 codons...")
        
        # Initialize counts for all possible codons
        all_codons = []
        bases = ['A', 'T', 'G', 'C']
        for base1 in bases:
            for base2 in bases:
                for base3 in bases:
                    all_codons.append(base1 + base2 + base3)
        
        self.all_codon_counts = {codon: 0 for codon in all_codons}
        
        # Get CDS features from annotation
        cds_features = self.annotation_parser.stats.get('cds_features', [])
        
        if not cds_features:
            print("⚠️  No CDS features found for all-codon analysis")
            return
        
        # Get genome sequences
        sequences = {}
        for seq_record in self.genome_parser.sequences:
            sequences[seq_record.id] = str(seq_record.seq).upper()
        
        total_codons_found = 0
        cds_analyzed = 0
        
        for cds in cds_features:
            seqid = cds['seqid']
            start = cds['start'] - 1  # Convert to 0-based indexing
            end = cds['end']
            strand = cds['strand']
            
            if seqid not in sequences:
                continue
            
            # Extract CDS sequence
            cds_seq = sequences[seqid][start:end]
            
            # Handle reverse strand
            if strand == '-':
                cds_seq = self._reverse_complement(cds_seq)
            
            # Count all codons in this CDS (in triplets)
            for i in range(0, len(cds_seq) - 2, 3):
                codon = cds_seq[i:i+3]
                if len(codon) == 3 and all(base in 'ATGC' for base in codon):
                    if codon in self.all_codon_counts:
                        self.all_codon_counts[codon] += 1
                        total_codons_found += 1
            
            cds_analyzed += 1
        
        print(f"   📊 Analyzed {cds_analyzed} CDS sequences")
        print(f"   📊 Found {total_codons_found:,} total codons")
        print("✅ All-codon analysis complete")

class ExonAnalyzer:
    """Analyzes Coding Sequence Exons - CDS regions from gene annotations"""
    
    def __init__(self, genome_parser: GenomeParser, annotation_parser: AnnotationParser):
        self.genome_parser = genome_parser
        self.annotation_parser = annotation_parser
        self.exon_stats = {}
    
    def _extract_gene_info(self, attributes: str) -> Dict[str, str]:
        """Extract gene name and product from GFF attributes"""
        info = {'gene': 'Unknown', 'product': 'Unknown'}
        
        # Parse attributes string (format: key=value;key=value;...)
        for attr in attributes.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                key = key.strip().lower()
                value = value.strip()
                
                # Clean up values (remove quotes and URL encoding)
                if value.startswith('"') and value.endswith('"'):
                    value = value[1:-1]
                value = value.replace('%2C', ',').replace('%3B', ';')
                
                # Extract gene names
                if key in ['gene', 'gene_name', 'name', 'locus_tag']:
                    if info['gene'] == 'Unknown' or key == 'gene':  # Prefer 'gene' over others
                        info['gene'] = value
                # Extract product/function information
                elif key in ['product', 'note', 'function', 'description']:
                    if info['product'] == 'Unknown' or key == 'product':  # Prefer 'product' over others
                        info['product'] = value
        
        # If gene is still unknown but we have an ID-like product, use a shorter version
        if info['gene'] == 'Unknown' and info['product'] != 'Unknown':
            product_words = info['product'].split()
            if len(product_words) > 0 and any(c.isdigit() for c in product_words[0]):
                info['gene'] = product_words[0]  # Use first word if it contains numbers (likely an ID)
        
        return info

    def analyze_exons(self) -> Dict:
        """Analyze Exon (CDS) statistics based on real CDS sequences"""
        print("🔬 Analyzing Exons (Coding Sequence regions)...")
        
        # Get CDS features from annotation
        cds_features = self.annotation_parser.stats.get('cds_features', [])
        
        if not cds_features:
            print("⚠️  No CDS features found in annotation file")
            print("❌ Cannot perform Exon analysis without CDS annotations")
            # Return empty stats instead of mock data
            self.exon_stats = {
                'total_exons': 0,
                'mean_length': 0,
                'median_length': 0,
                'longest_exon': 0,
                'shortest_exon': 0
            }
            return self.exon_stats
        
        # Calculate Exon lengths from CDS features and track longest/shortest
        exon_data = []
        for cds in cds_features:
            length = cds['end'] - cds['start'] + 1
            gene_info = self._extract_gene_info(cds.get('attributes', ''))
            exon_data.append({
                'length': length,
                'gene': gene_info['gene'],
                'product': gene_info['product'],
                'start': cds['start'],
                'end': cds['end']
            })
        
        if not exon_data:
            print("⚠️  No valid Exon lengths found")
            print("❌ Cannot calculate Exon statistics")
            # Return empty stats instead of mock data
            self.exon_stats = {
                'total_exons': 0,
                'mean_length': 0,
                'median_length': 0,
                'longest_exon': 0,
                'shortest_exon': 0
            }
            return self.exon_stats
        
        # Find longest and shortest Exons
        longest_exon = max(exon_data, key=lambda x: x['length'])
        shortest_exon = min(exon_data, key=lambda x: x['length'])
        
        # Calculate statistics
        import statistics
        exon_lengths = [exon['length'] for exon in exon_data]
        
        self.exon_stats = {
            'total_exons': len(exon_lengths),
            'mean_length': round(statistics.mean(exon_lengths), 1),
            'median_length': statistics.median(exon_lengths),
            'longest_exon': max(exon_lengths),
            'shortest_exon': min(exon_lengths),
            'longest_exon_info': longest_exon,
            'shortest_exon_info': shortest_exon
        }
        
        print(f"   🧬 Analyzed {len(exon_lengths)} Exons")
        print(f"   📏 Length range: {min(exon_lengths)} - {max(exon_lengths)} bp")
        print("✅ Exon analysis complete")
        return self.exon_stats
    
    def find_orfs_in_sequence(self, sequence: str, min_length: int = 100) -> List[Dict]:
        """Find ORFs in a DNA sequence (alternative method for annotation-free analysis)"""
        orfs = []
        start_codons = ['ATG', 'GTG', 'TTG']
        stop_codons = ['TAA', 'TAG', 'TGA']
        
        # Check all three reading frames
        for frame in range(3):
            i = frame
            while i < len(sequence) - 2:
                # Look for start codon
                codon = sequence[i:i+3]
                if codon in start_codons:
                    start_pos = i
                    # Look for stop codon
                    for j in range(i + 3, len(sequence) - 2, 3):
                        stop_codon = sequence[j:j+3]
                        if stop_codon in stop_codons:
                            orf_length = j + 3 - start_pos
                            if orf_length >= min_length:
                                orfs.append({
                                    'start': start_pos,
                                    'end': j + 3,
                                    'length': orf_length,
                                    'frame': frame + 1,
                                    'start_codon': codon,
                                    'stop_codon': stop_codon
                                })
                            i = j + 3
                            break
                    else:
                        # No stop codon found, move to next position
                        i += 3
                else:
                    i += 3
        
        return orfs

class ReportGenerator:
    """Generates output reports"""
    
    def __init__(self):
        self.results = {}
    
    def generate_summary(self, genome_stats: Dict, codon_stats: Dict, exon_stats: Dict) -> str:
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

📏 Exon Statistics:
- Total Exons: {exon_stats['total_exons']:,}
- Mean length: {exon_stats['mean_length']:.1f} bp
- Median length: {exon_stats['median_length']} bp"""
        
        # Add detailed info for longest and shortest Exons if available
        if 'longest_exon_info' in exon_stats and exon_stats['longest_exon_info']:
            longest_info = exon_stats['longest_exon_info']
            gene_name = longest_info.get('gene', 'Unknown')
            product = longest_info.get('product', 'Unknown')
            summary += f"\n- Longest Exon: {exon_stats['longest_exon']:,} bp"
            summary += f"\n  📍 Gene: {gene_name}"
            if product != 'Unknown' and product != gene_name:
                summary += f"\n  🔬 Product: {product}"
        else:
            summary += f"\n- Longest Exon: {exon_stats['longest_exon']:,} bp"
        
        if 'shortest_exon_info' in exon_stats and exon_stats['shortest_exon_info']:
            shortest_info = exon_stats['shortest_exon_info']
            gene_name = shortest_info.get('gene', 'Unknown')
            product = shortest_info.get('product', 'Unknown')
            summary += f"\n- Shortest Exon: {exon_stats['shortest_exon']} bp"
            summary += f"\n  📍 Gene: {gene_name}"
            if product != 'Unknown' and product != gene_name:
                summary += f"\n  🔬 Product: {product}"
        else:
            summary += f"\n- Shortest Exon: {exon_stats['shortest_exon']} bp"
        
        summary += f"""

🍽️ Nom nom nom! Analysis complete.
"""
        return summary
    
    def save_csv(self, filename: str, codon_stats: Dict, exon_stats: Dict):
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
            
            # Add Exon stats
            for metric, value in exon_stats.items():
                writer.writerow({
                    'metric_type': 'exon_stats',
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
  python genomnomnom.py --genome genome.fasta --annotation genes.gff --detailed-codons
  
  # Search and download from NCBI
  python genomnomnom.py --species "Escherichia coli" --email your.email@example.com
  python genomnomnom.py --species "Homo sapiens" --email your.email@example.com --output human_results.csv
  python genomnomnom.py --species "Escherichia coli" --email your.email@example.com --detailed-codons
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
    
    parser.add_argument(
        "--detailed-codons", "-d",
        action="store_true",
        help="Show detailed codon usage table with all 64 codons"
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
        
        # Analyze Exons
        exon_analyzer = ExonAnalyzer(genome_parser, annotation_parser)
        exon_stats = exon_analyzer.analyze_exons()
        
        # Generate report
        report_generator = ReportGenerator()
        summary = report_generator.generate_summary(genome_stats, codon_stats, exon_stats)
        
        print(summary)
        
        # Show detailed codon table if requested
        if args.detailed_codons:
            codon_table = codon_analyzer.generate_codon_table_report()
            print(codon_table)
            
            # Also show classic square table
            classic_table = codon_analyzer.generate_classic_codon_table()
            print(classic_table)
        
        # Save CSV if requested
        if args.output:
            report_generator.save_csv(args.output, codon_stats, exon_stats)
        
        print("🎉 GenomNomNom finished successfully!")
        
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
