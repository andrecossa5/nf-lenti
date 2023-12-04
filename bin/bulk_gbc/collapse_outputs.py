#!/usr/bin/python

# Create summary output script 

########################################################################

# Parsing CLI args 

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
	prog='collapse_outputs',
	description=
	"""
	Script for aggregating outputs nf-perturbseq pipeline, bulk.
	"""
)

# Add arguments
 
# Input
my_parser.add_argument(
	'-i', 
	'--input', 
	type=str,
	default=None,
	help='The path to results/perturb_bulk folder. Default: None.'
)
 
# Ouput
my_parser.add_argument(
	'-o', 
	'--output', 
	type=str,
	default=os.getcwd(),
	help='''
		The path to the output folder of choice. Default: os.getcwd().
		'''
)
 
 
# Parse arguments
args = my_parser.parse_args()
path_input = args.input
path_output = args.output

##

# Import code
import numpy as np
import pandas as pd

##

##########################################################


# Main
def main():

	# Create summary dir
	os.mkdir('summary')
	os.chdir('summary')

	# Read all sets
	d = {}
	for root, _, files in os.walk(path_input):
	    for file in files:
	        if file == 'GBC_counts_corrected.csv':
	            sample = os.path.basename(root)
	            d[sample] = (
					pd.read_csv(os.path.join(root, file), index_col=0)
					.assign(sample=sample)
				)

	# Final reference
	pd.concat(d.values()).to_csv('bulk_GBC_reference.csv')

	# Common clones
	n = len(d)
	C = np.zeros((n,n))
	for i,x in enumerate(d.keys()):
		for j,y in enumerate(d.keys()):
			if i>=j:
				C[i,j] = C[j,i] = len(set(d[x].index)&set(d[y].index))		
	pd.DataFrame(C, index=d.keys(), columns=d.keys()).to_csv('common.csv')
	
	##

##########################################################

# Run 
if __name__ == '__main__':
	main()
