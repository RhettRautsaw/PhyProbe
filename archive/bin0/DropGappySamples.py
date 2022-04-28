#!/usr/bin/env python3

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import multiprocessing as mp
import glob
import pandas as pd
import numpy as np
from tqdm import tqdm
from p_tqdm import p_map
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

Blah blah blah

""")

parser.add_argument("-f","--folder",
					type=str,
					default="05_trimal",
					help="Folder of alignments")
parser.add_argument("-o","--output",
					type=str,
					default="DropGappySeqs",
					help="Folder for alignments without gappy samples")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate target/gene names from sample names (default: %(default)s)")
parser.add_argument("-g","--gaps",
					type=int,
					default=25,
					help="Percent of sequence allowed to be gaps (default: %(default)s)")
parser.add_argument("-c","--cpu",
					type=int,
					default=mp.cpu_count(),
					help="Number of threads to be used in each step. (default: %(default)s)")
args=parser.parse_args()

########################################
############### FUNCTIONS ##############
########################################

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc: # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else: raise

def MissingData(file):
	locus = file.split(".")[0].split("/")[-1]
	tmp_seq = list(SeqIO.parse(file,'fasta'))
	clean_seqs = []
	data = dict()
	for seq in tmp_seq:
		sample = seq.name.split("|")[0]
		c = seq.seq.count("-")
		l = len(seq.seq)
		completeness = 1-(c/l)
		if completeness > (1-gaps):
			clean_seqs.append(seq)
			data[sample] = completeness
	
	handle=open(os.path.join(output, locus) + ".fasta",'w')
	writer = FastaIO.FastaWriter(handle)
	writer.write_file(clean_seqs)
	handle.close()

########################################
################# SETUP ################
########################################

if args.folder==None:
	print("No input folder specified. Use -f/--folder flag to specify")
	quit()
else:
	folder = os.path.abspath(args.folder)
	files = glob.glob(folder + "/*.fasta")

if args.output==None:
	print("No output location specified. Use -o/--output flag to specify")
	quit()
else:
	output = args.output
	mkdir_p(output)

delim=args.delim
gaps = args.gaps/100
cpus = args.cpu

print("""

Blah blah blah

""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting...")
print("\tInput folder -> "+ args.folder)
print("\tOutput folder -> "+ output)
print("\tDelimiter -> " + delim)
print("\tPercent Gaps Allowed in Sequence -> " + str(args.gaps) + "%")
print("\tNumber of CPU -> "+ str(cpus))

########################################
################# CODE #################
########################################

p_map(MissingData, files, num_cpus=cpus)
