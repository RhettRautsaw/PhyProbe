#!/usr/bin/env python3

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import gzip
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

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

For duplicate genes/isoforms, this script will figure out which is a the longest and only keep that one.

Default delimiter ( -d | ) and position ( -p 0 ) arguments assume this format:
>gene_name|sample_name|any_other_information

For Trinity isoforms:
	-d _i -p 0

For rnaSPADES isoforms:
	-d _ -p 7  or -2

""")

########################################
############### ARGUMENTS ##############
########################################

parser.add_argument("-i","--input",
					type=str,
					help="File with sequences to check agains target database (default: %(default)s)")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate target/gene names from sample names (default: %(default)s)")
parser.add_argument("-p","--position",
					type=int,
					default=0,
					help="After splitting by delimiter, what 0-indexed position is used to identify the gene/isoform (default: %(default)s)")
parser.add_argument("--auto",
					action="store_true",
					default=False,
					help="Add this flag to auto-detect Trinity and rnaSpades assemblies for isoform parsing (default: %(default)s)")
parser.add_argument("-o","--output",
					type=str,
					default='.',
					help="Folder in which to export reduced fasta with only longest forms (default: %(default)s)")
args=parser.parse_args()

########################################
################# SETUP ################
########################################

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc: # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else: raise

input=os.path.abspath(args.input)
outfile=input.split("/")[-1]
outfile=outfile.split(".")[:-1]
outfile='.'.join(outfile)+".long.fasta"

mkdir_p(args.output)
outfile = os.path.join(os.path.abspath(args.output), outfile)

delim = args.delim
pos = args.position

########################################
################# START ################
########################################

fasta = list(SeqIO.parse(input,'fasta'))

if args.auto:
	if 'NODE' in fasta[0].name:
		print("rnaSpades assembly detected")
		delim='_'
		pos=7
	elif 'TRINITY' in fasta[0].name:
		print("Trinity assembly detected")
		delim='_i'
		pos=0
	else:
		print("Trinity or rnaSpades assembly not detected")

names=[]
genes=[]
lengths=[]

for seq in fasta:
	gene = seq.name.split(delim)[pos]
	names.append(seq.name)
	genes.append(gene)
	lengths.append(len(seq.seq.ungap("-")))

df = {'transcript': names, 'gene': genes, 'length': lengths}
df = pd.DataFrame(df) 

df2 = df >> group_by(X.gene) >> arrange(X.length, ascending=False) >> head(1)
keep = list(df2.transcript)

long_locus = []
for seq in fasta:
	if seq.name in keep:
		long_locus.append(seq)


handle=open(outfile, "w")
writer = FastaIO.FastaWriter(handle)
writer.write_file(long_locus)
handle.close()
