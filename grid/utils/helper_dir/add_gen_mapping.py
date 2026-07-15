#!/usr/bin/env python3

__author__ = 'caterer-z-t'

import pandas as pd
import numpy as np
import argparse
import sys

def main(args):
    # Read the Eagle genetic map
    genetic_map = pd.read_csv(
        args.genetic_map,
        sep='\s+',
        compression='gzip',
        comment='#'
    )
    
    # Read MAP file
    map_file = pd.read_csv(
        args.map, 
        sep='\t', 
        names=['chr', 'snp', 'cm', 'bp']
    )
    
    # Interpolate genetic positions for MAP file variants
    map_file['cm'] = np.interp(
        map_file['bp'],
        genetic_map['position'],
        genetic_map['Genetic_Map(cM)']
    )
    
    # Save updated MAP file
    output_file = f"{args.out}.map"
    map_file.to_csv(
        output_file, 
        sep='\t', 
        header=False, 
        index=False
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add genetic distances to PLINK MAP file using Eagle genetic map')
    
    parser.add_argument('--map', help='PLINK MAP file', required=True)
    parser.add_argument('--genetic-map', help='Eagle genetic map file (gzipped)', required=True)
    parser.add_argument('--out', help='prefix for output', required=True)
    
    args = parser.parse_args()
    main(args)
    sys.exit()