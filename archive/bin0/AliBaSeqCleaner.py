#!/usr/bin/env python3

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import multiprocessing as mp
import glob
from dfply import *
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

AliBaSeq Cleaner

Results will have the format:
>locus_name|sample_name

Delimiter | can be changed with the -d flag.

""")

parser.add_argument("-f","--folder",
					type=str,
					default="alibaseq_out",
					help="AliBaSeq output folder")
parser.add_argument("-n","--name",
					type=str,
					help="Name to apply to sequence")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate locus names from sample names (default: %(default)s)")
parser.add_argument("--percent",
					type=int,
					default=25,
					help="Percent missing data (gaps) allowed in a locus. (default: %(default)s)")
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

def Cleaner(file):
	locus = file.split("/")[-1].split(".fas")[0]
	tmp_seq = list(SeqIO.parse(file,'fasta'))
	for seq in tmp_seq:
		c = seq.seq.count("N")
		l = len(seq.seq)
		if c/l < percent:
			seq.seq=Seq(str(seq.seq).replace("N","-"))
			seq.name = seq.id = seq.description = locus + delim + name
			output_handle = open(output + "/genes/" + locus + ".fasta",'w')
			SeqIO.write(seq, output_handle, "fasta")
			output_handle.close()
			sp.call('cat '+ output + '/genes/' + locus + '.fasta >> ' + os.path.join(output,name) + '.targets.fasta', shell=True)

########################################
################# SETUP ################
########################################

folder = os.path.abspath(args.folder)
output = '/'.join(folder.split("/")[:-1])
output = os.path.join(output, "alibaseq_clean")

if args.name==None:
	print("No sample name specified. Use -n/--name flag to specify input")
	quit()
else:
	name = args.name

delim = args.delim
percent = args.percent/100
cpus = args.cpu

print("""

Clean AliBaSeq Results

Results will have the format:
>locus_name|sample_name

Delimiter | can be changed with the -d flag.

""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting AliBaSeq Cleaner...")
print("\tAliBaSeq output folder -> "+ args.folder)
print("\tSample name -> "+ args.name)
print("\tDelimiter -> " + delim)
print("\tPercent Missing Data -> " + str(args.percent))
print("\tNumber of CPU -> "+ str(cpus))

########################################
################# CODE #################
########################################

files = glob.glob(os.path.join(folder,"*.fas"))

mkdir_p(output)
mkdir_p(output + "/genes")
results=p_map(Cleaner, files, num_cpus=cpus)
