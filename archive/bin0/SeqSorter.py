#!/usr/bin/env python3

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import gzip
from tqdm import tqdm
from p_tqdm import p_map
from dfply import *
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.SeqIO import FastaIO
except:
	print("Error: biopython module is not properly installed.")
	quit()

try:
	import numpy as np
except:
	print("Error: numpy is not properly installed.")
	quit()

try:
	import pandas as pd
except:
	print("Error: pandas is not properly installed.")
	quit()

try:
	import glob
except:
	print("Error: glob is not properly installed.")
	quit()

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

Sort sequences into separate files by sample or locus

Assumes targets have the format:
>target_name|sample_name

Delimiter | can be changed with the -d flag.
Position can be changed with the -p flag.

Default is -d "|" and -p 1

""")

########################################
############### ARGUMENTS ##############
########################################

parser.add_argument("-i","--input",
					type=str,
					help="File with sequences to sort by sample or gene name (default: %(default)s)")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate target/gene names from sample names (default: %(default)s)")
parser.add_argument("-p","--position",
					type=int,
					default=1,
					help="0-indexed position to determine sorting by sample or locus. (default: %(default)s)")
parser.add_argument("-c","--cpu",
					type=int,
					default=8,
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

def SeqSorter(seq):
	sample=seq.name.split(delim)[pos]
	
	output_handle = open(folder + "/SeqSorter/" + sample + '.fasta','a')
	SeqIO.write(seq, output_handle, "fasta")
	output_handle.close()

########################################
################# SETUP ################
########################################

input=os.path.abspath(args.input)
folder=input.split("/")[:-1]
folder='/'.join(folder)
delim = args.delim
pos = args.position
cpus = args.cpu

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting SeqSorter...")
print("\tInput -> "+ input)
print("\tDelim -> "+ delim)
print("\tThreads -> "+ str(cpus))

########################################
################# CODE #################
########################################

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Reading ' + args.input)
sequences = list(SeqIO.parse(input,"fasta"))

d={}
print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Parsing sequences into samples :::')
for seq in tqdm(sequences):
	sample=seq.name.split(delim)[pos]
	if sample not in d:
		d.setdefault(sample, [])
		d[sample].append(seq)
	else:
		d[sample].append(seq)

mkdir_p(folder + '/SeqSorter2')
for locus in tqdm(d):
	tmp_list=d[locus]
	output_handle = open(folder + "/SeqSorter2/" + locus + '.fasta','w')
	SeqIO.write(tmp_list, output_handle, "fasta")
	output_handle.close()

#results = p_map(SeqSorter, sequences, num_cpus=cpus)
