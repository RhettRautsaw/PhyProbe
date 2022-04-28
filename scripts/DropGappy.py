#!/usr/bin/env python

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
from tqdm import tqdm
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.SeqIO import FastaIO
except:
	print("Error: biopython module is not properly installed.")
	quit()

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

Get rid of those pesky samples/sequences that don't have good coverage of the alignment.

""")

parser.add_argument("-i","--input",
					type=str,
					help="Input alignment")
parser.add_argument("-o","--output",
					type=str,
					help="Output alignment")
parser.add_argument("-c","--coverage",
					type=int,
					default=50,
					help="Percent of alignment coverage required (i.e., percent of sample that is not gaps; default: %(default)s)")
parser.add_argument("-ug","--ungap",
					action="store_true",
					default=False,
					help="Remove all gaps from alignment so that you can realign without gappy sequences? (default: %(default)s)")
parser.add_argument("-w","--wrap",
					action="store_true",
					default=False,
					help="Default is to produce linear fasta, but this flag will wrap your fastas for you.")
args=parser.parse_args()

########################################
################# SETUP ################
########################################

if args.input==None:
	print("No input alignment specified. Use -i/--input flag to specify")
	quit()
else:
	input = os.path.abspath(args.input)
	input_basename = os.path.basename(input).split('.fa')[0]

if args.output==None:
	output=input_basename+"_dropgappy.fasta"
else:
	output = args.output


req_cov = args.coverage/100

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting...")
print("\tInput alignment -> "+ input)
print("\tOutput alignment -> "+ output)
print("\tRequired Coverage (non-gaps) -> " + str(args.coverage) + "%")
print("\tRemove Gaps from Output -> " + str(args.ungap))
print("\tWrap Fasta Output -> " + str(args.wrap))

########################################
################# CODE #################
########################################

locus = list(SeqIO.parse(input,'fasta'))
clean_seqs = []
for seq in tqdm(locus):
	c = seq.seq.count("-")
	l = len(seq.seq)
	coverage=(1-(c/l))
	if coverage > req_cov:
		if args.ungap:
			seq.seq==seq.seq.ungap("-")
		clean_seqs.append(seq)

print("Writing " + str(len(clean_seqs)) + " of " + str(len(locus)) + " total sequences")

## Write results to file
if args.wrap==False:
	SeqIO.write(clean_seqs, output, "fasta-2line")
if args.wrap==True:
	SeqIO.write(clean_seqs, output, "fasta")
