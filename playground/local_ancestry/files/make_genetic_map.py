#!/usr/bin/env python
'''
creates a genetic map of sites from a vcf using sites in a pre-made genetic map
'''
__author__ = 'jashortt'


from numpy import *
from scipy.interpolate import interp1d
import argparse
import sys
import gzip
import random				# might not be used
import os					# might not be used

def printDict(hash):
	for x in hash:
		print (x, hash[x])
	assert 1 == 0

def getDistsFromMap (mapfile):
	print ("Getting distances from file " + mapfile)
	dists = {}
	with gzip.open(mapfile, 'rt') as infp:
		header = infp.readline()
		for line in infp:
			chrom, pos, rate, dist = tuple(line.strip().split())
			if not chrom in dists:
				dists[chrom] = {}	
			chrom_dists = dists[chrom]
			if not 'dist' in chrom_dists:
				chrom_dists['dist'] = []
			if not 'pos' in chrom_dists:
				chrom_dists['pos'] = []
			chrom_dists['dist'].append (float(dist))
			chrom_dists['pos'].append (int (pos))
	print ('Finished getting distances from file ' + mapfile)
	return (dists)
	
def getSitesFromVcf (vcf, mychrom):
	sites = []
	print ('Getting sites of interest from file ' + vcf)
	with gzip.open(vcf, 'rt') as infp:
		for line in infp:
			if line.startswith("#"):
				pass
			else:
				info = line.strip().split() 
				chrom, pos = tuple(info[:2])
				if mychrom != chrom:
					print ('Error:\n' + chrom + ' is not the same chromosome as ' + mychrom)
					assert 1 == 0
				sites.append( int(pos) )
	return sites
        
def printNewMap (file, chrom, sites, bpcms):
	print ("Printing site map distances to " + file)
	with open (file, 'w') as outfp:
		for pos, dist in zip(sites,bpcms):
			outfp.write( "%s\t%s\t%8.5f\n" % (str(chrom), str(pos), dist) )
		
def main (args):
	dists = getDistsFromMap (args.map)
	sites = getSitesFromVcf (args.vcf, args.chrom)
	# create interpolation function
	f = interp1d(dists[args.chrom]['pos'],dists[args.chrom]['dist'],fill_value='extrapolate')
	bpcms = f(sites)
	printNewMap("%s.gnomix.map" % (args.out), args.chrom, sites, bpcms)
	
		
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--vcf', help='zipped .vcf file', required=True)
    parser.add_argument('--map', help='zipped map .txt from Eagle (chrom\tpos\trate\tdist)', required=True)
    parser.add_argument('--out', help='prefix for output', required=True)
    parser.add_argument('--chrom', help='chromosome of interest', required=True)

    args = parser.parse_args()
    main(args)
    sys.exit()