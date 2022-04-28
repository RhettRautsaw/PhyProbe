#!/usr/bin/env python

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import glob
from tqdm import tqdm
import pyfastx
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.SeqIO import FastaIO
except:
	print("Error: biopython module is not properly installed.")
	quit()

import multiprocessing as mp
#from dfply import *
#from p_tqdm import p_map

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
					help="AliBaSeq results folder")
parser.add_argument("-r","--reference",
					type=str,
					default="",
					help="Reference sequence fasta")
parser.add_argument("-n","--name",
					type=str,
					help="Name to apply to sequence")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate locus names from sample names (default: %(default)s)")
parser.add_argument("--percent",
					type=int,
					default=50,
					help="Percent coverage of reference locus required. (default: %(default)s)")
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

#def Cleaner(file):
#	locus = file.split("/")[-1].split(".fas")[0]
#	tmp_seq = list(SeqIO.parse(file,'fasta'))
#	reference = pyfastx.Fasta(reference_name)
#	tl=len(reference[locus])
#	for seq in tmp_seq:
#		l = len(seq.seq)-seq.seq.count("N")
#		if l/tl > percent:
#			seq.seq=Seq(str(seq.seq).replace("N","-"))
#			seq.name = seq.id = seq.description = locus + delim + name
#			output_handle = open(output + "/genes/" + locus + ".fasta",'w')
#			SeqIO.write(seq, output_handle, "fasta")
#			output_handle.close()
#			sp.call('cat '+ output + '/genes/' + locus + '.fasta >> ' + os.path.join(output,name) + '_targets.fasta', shell=True)

########################################
################# SETUP ################
########################################

folder = os.path.abspath(args.folder)
reference_name = os.path.abspath(args.reference)
reference = pyfastx.Fasta(reference_name)

ref_dic={}
for i in reference:
	if i.name not in ref_dic.keys():
		ref_dic[i.name]=len(i)
	else:
		if ref_dic[i.name] < len(i):
			ref_dic[i.name]=len(i)

if args.name==None:
	print("No sample name specified. Use -n/--name flag to specify input")
	quit()
else:
	name = args.name

delim = args.delim
percent = args.percent/100
cpus = args.cpu

output = "alibaseq_clean"
output_name=os.path.join(output,name) + '_targets.fasta'

print("""

Clean AliBaSeq Results

Results will have the format:
>locus_name|sample_name

Delimiter | can be changed with the -d flag.

""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting AliBaSeq Cleaner...")
print("\tAliBaSeq Results -> "+ folder)
print("\tReference Fasta -> "+ reference_name)
print("\tSample name -> "+ args.name)
print("\tDelimiter -> " + delim)
print("\tReference Locus Coverage -> " + str(args.percent) + "%")
print("\tNumber of CPU -> "+ str(cpus))
print("\tOutput -> "+ output_name)

########################################
################# CODE #################
########################################

files = glob.glob(os.path.join(folder,"*.fas"))

mkdir_p(output)
#mkdir_p(output + "/genes")
#results=p_map(Cleaner, files, num_cpus=cpus)

output_handle = open(output_name,'a')

for file in tqdm(files):
	locus = file.split("/")[-1].split(".fas")[0]
	tmp_seq = list(SeqIO.parse(file,'fasta'))
	tl=ref_dic[locus]
	for seq in tmp_seq:
		l = len(seq.seq)-seq.seq.count("N")
		if l/tl > percent:
			seq.seq=Seq(str(seq.seq).replace("N","-"))
			seq.name = seq.id = seq.description = locus + delim + name
			SeqIO.write(seq, output_handle, "fasta")

output_handle.close()

new_seqs = list(SeqIO.parse(output_name,'fasta'))

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: " + str(len(files)) + " sequences reduced to " + str(len(new_seqs)) + " :::\n")
