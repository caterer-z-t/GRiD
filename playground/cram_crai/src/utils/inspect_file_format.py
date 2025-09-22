#!/usr/bin/env python3
"""
Inspect the format of the normalized mosdepth file to understand its structure
"""

import gzip
import sys

def inspect_file(filename):
    """Inspect the first few lines of the file to understand its format"""
    
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'rt')
    else:
        f = open(filename, 'r')
    
    print(f"Inspecting file: {filename}")
    print("=" * 50)
    
    try:
        # Read first 10 lines
        for i in range(10):
            line = f.readline()
            if not line:
                break
            line = line.strip()
            print(f"Line {i+1}: {line}")
            
            # Try to analyze the structure
            if i == 0:
                parts = line.split()
                print(f"  -> {len(parts)} parts: {parts[:10]}...")  # Show first 10 parts
                
        print("\n" + "=" * 50)
        print("Analysis:")
        
        # Go back to start and analyze structure
        f.seek(0)
        first_line = f.readline().strip()
        parts = first_line.split()
        
        print(f"First line has {len(parts)} parts")
        
        # Check if first two parts are integers (N and R)
        try:
            n_val = int(parts[0])
            r_val = int(parts[1])
            print(f"Looks like standard format: N={n_val}, R={r_val}")
        except:
            print("First line doesn't start with two integers")
            print("This might be a header line or different format")
            
            # Check if it starts with sample ID
            if parts[0].startswith('HG') or parts[0].startswith('NA') or len(parts[0]) > 5:
                print("First part looks like a sample ID - this might be data, not header")
                print("Expected format: N R on first line, but got data instead")
    
    except Exception as e:
        print(f"Error reading file: {e}")
    
    finally:
        f.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python inspect_file_format.py <input_file>")
        sys.exit(1)
    
    inspect_file(sys.argv[1])