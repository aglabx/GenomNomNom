#!/usr/bin/env python3
"""
Test script for NCBI search functionality
"""

import sys
sys.path.insert(0, '.')
from genomnomnom import NCBISearcher

def test_ncbi_search():
    """Test NCBI search without user interaction"""
    print("🧪 Testing NCBI search functionality...")
    
    ncbi = NCBISearcher(email="test@example.com")
    
    # Test search
    species = "Mycoplasma genitalium"
    print(f"🔍 Searching for '{species}'...")
    
    assemblies = ncbi.search_genomes(species, max_results=5)
    
    if assemblies:
        print(f"✅ Found {len(assemblies)} assemblies")
        print("\n📋 Available assemblies:")
        for i, assembly in enumerate(assemblies, 1):
            size_mb = int(assembly['size']) / 1_000_000 if assembly['size'] else 0
            print(f"{i:2d}. {assembly['organism']}")
            print(f"    Accession: {assembly['accession']}")
            print(f"    Name: {assembly['name']}")
            print(f"    Level: {assembly['level']}")
            print(f"    Size: {size_mb:.1f} MB")
            print()
        
        # Test with the first assembly (without download)
        print("🔬 Testing download URLs construction...")
        first_assembly = assemblies[0]
        if first_assembly.get('ftp_path'):
            accession = first_assembly['accession']
            ftp_base = first_assembly['ftp_path'].replace('ftp://', 'https://')
            base_name = ftp_base.split('/')[-1]
            fasta_url = f"{ftp_base}/{base_name}_genomic.fna.gz"
            gff_url = f"{ftp_base}/{base_name}_genomic.gff.gz"
            
            print(f"   FASTA URL: {fasta_url}")
            print(f"   GFF URL: {gff_url}")
            
            # Test if URLs are reachable (HEAD request)
            import requests
            try:
                response = requests.head(fasta_url, timeout=10)
                print(f"   FASTA file exists: {response.status_code == 200}")
            except:
                print("   FASTA file check failed")
                
            try:
                response = requests.head(gff_url, timeout=10)
                print(f"   GFF file exists: {response.status_code == 200}")
            except:
                print("   GFF file check failed")
        
        print("✅ NCBI search test completed successfully!")
    else:
        print("❌ No assemblies found")

if __name__ == "__main__":
    test_ncbi_search()
