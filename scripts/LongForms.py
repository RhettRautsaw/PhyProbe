#!/usr/bin/env python

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import gzip
#from dfply import *
import numpy as np
import pandas as pd
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.SeqIO import FastaIO
except:
	print("Error: biopython module is not properly installed.")
	quit()


parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

Identify and Extract the longest isoform from an RNA assembly.

For Trinity isoforms:
	-d _i -p 0

For rnaSPADES isoforms:
	-d _  -p 6

""")

########################################
############### ARGUMENTS ##############
########################################

parser.add_argument("-i","--input",
					type=str,
					default=None,
					help="Assembly fasta (default: %(default)s)")
parser.add_argument("--manual",
					action="store_true",
					default=False,
					help="LongForms.py will auto-detect Trinity and rnaSpades assemblies by default, but you can add this flag to manually input delimiters and position arguments. (default: %(default)s)")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Optional manual delimiter used to isolate gene names in fasta header (default: %(default)s)")
parser.add_argument("-p","--position",
					type=int,
					default=0,
					help="Optional manual position argument for 0-indexed position of gene names in fasta header (default: %(default)s)")
parser.add_argument("-o","--output",
					type=str,
					default=None,
					help="Name of file to output (default: %(default)s)")
args=parser.parse_args()

if args.input==None:
	print("No input fasta provided. Use -i/--input flag to specify input.")
	quit()
else:
	input = os.path.abspath(args.input)

if args.output==None:
	print("No output name given. Use -o/--output flag to provide file name for output.")
	quit()
else:
	output = os.path.abspath(args.output)


########################################
################# START ################
########################################

fasta = list(SeqIO.parse(input,'fasta'))

if args.manual:
	delim = args.delim
	pos = args.position
else:
	if 'NODE' in fasta[0].name:
		print("rnaSpades assembly detected")
		delim='_'
		pos=6
	elif 'TRINITY' in fasta[0].name:
		print("Trinity assembly detected")
		delim='_i'
		pos=0
	else:
		print("Neither trinity or rnaSpades assembly detected. Use --manual")
		quit()

names=[]
genes=[]
lengths=[]

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Building Database :::')
for seq in fasta:
	gene = seq.name.split(delim)[pos]
	names.append(seq.name)
	genes.append(gene)
	lengths.append(len(seq.seq.ungap("-")))

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: ID LongForm :::')
df = {'transcript': names, 'gene': genes, 'length': lengths}
df = pd.DataFrame(df) 
#df2 = df >> group_by(X.gene) >> arrange(X.length, ascending=False) >> head(1)
df2 = df.sort_values("length", ascending=False).groupby("gene").head(1)
keep = list(df2.transcript)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Extract LongForm :::')
long_locus = []
for seq in fasta:
	if seq.name in keep:
		long_locus.append(seq)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Writing Results :::')
handle=open(args.output, "w")
writer = FastaIO.FastaWriter(handle)
writer.write_file(long_locus)
handle.close()
