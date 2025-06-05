#!/usr/bin/env python3
"""
Test script to demonstrate GenomNomNom functionality
"""

import subprocess
import sys
from pathlib import Path

def test_genomnomnom():
    """Test the GenomNomNom tool with sample data"""
    
    # Get paths
    script_dir = Path(__file__).parent
    genome_file = script_dir / "test_data" / "genome_sample.fasta"
    annotation_file = script_dir / "test_data" / "annotation_sample.gff"
    output_file = script_dir / "test_results.csv"
    
    print("🧪 Testing GenomNomNom with sample data...")
    print("=" * 50)
    
    # Test basic functionality
    cmd = [
        sys.executable, "genomnomnom.py",
        "--genome", str(genome_file),
        "--annotation", str(annotation_file),
        "--output", str(output_file),
        "--verbose"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=script_dir)
        
        print("STDOUT:")
        print(result.stdout)
        
        if result.stderr:
            print("\nSTDERR:")
            print(result.stderr)
        
        if result.returncode == 0:
            print("✅ Test completed successfully!")
            
            # Check if output file was created
            if output_file.exists():
                print(f"📄 Output file created: {output_file}")
                with open(output_file, 'r') as f:
                    print("\nOutput file content preview:")
                    print(f.read()[:500] + "..." if len(f.read()) > 500 else f.read())
            else:
                print("⚠️  Output file was not created")
        else:
            print(f"❌ Test failed with return code: {result.returncode}")
            
    except Exception as e:
        print(f"❌ Error running test: {e}")

def test_help():
    """Test help functionality"""
    print("\n🆘 Testing help functionality...")
    print("=" * 30)
    
    script_dir = Path(__file__).parent
    cmd = [sys.executable, "genomnomnom.py", "--help"]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=script_dir)
        print(result.stdout)
    except Exception as e:
        print(f"❌ Error testing help: {e}")

if __name__ == "__main__":
    test_help()
    test_genomnomnom()
