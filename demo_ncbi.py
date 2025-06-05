#!/usr/bin/env python3
"""
Demo script for GenomNomNom NCBI functionality

This script demonstrates the new NCBI integration features without
requiring user interaction.
"""

import sys
import os
sys.path.insert(0, '.')

def demo_help():
    """Show the updated help message"""
    print("🎭 GenomNomNom NCBI Integration Demo")
    print("=" * 40)
    print()
    print("1. Let's see the updated help message:")
    print("-" * 40)
    
    os.system("python genomnomnom.py --help")
    print()

def demo_local_analysis():
    """Demo local file analysis"""
    print("2. Local file analysis (as before):")
    print("-" * 40)
    
    os.system("python genomnomnom.py --genome test_data/genome_sample.fasta --annotation test_data/annotation_sample.gff")
    print()

def demo_ncbi_test():
    """Demo NCBI functionality test"""
    print("3. NCBI search functionality test:")
    print("-" * 40)
    
    os.system("python test_ncbi.py")
    print()

def demo_summary():
    """Show summary of new features"""
    print("🎉 New Features Summary:")
    print("=" * 40)
    print("✅ NCBI Entrez API integration")
    print("✅ Species name search")
    print("✅ Interactive genome selection")
    print("✅ Automatic FASTA and GFF download")
    print("✅ Seamless integration with existing analysis")
    print()
    print("📚 Usage patterns:")
    print("  Local files:  python genomnomnom.py --genome file.fasta --annotation file.gff")
    print("  NCBI search:  python genomnomnom.py --species 'Species name' --email user@example.com")
    print()
    print("🧬 Happy genome munching!")

if __name__ == "__main__":
    print("🧬 GenomNomNom NCBI Integration Demo")
    print("🎯 This demo shows the new NCBI functionality")
    print()
    
    try:
        demo_help()
        demo_local_analysis()
        demo_ncbi_test()
        demo_summary()
    except KeyboardInterrupt:
        print("\n👋 Demo interrupted by user")
    except Exception as e:
        print(f"❌ Demo error: {e}")
