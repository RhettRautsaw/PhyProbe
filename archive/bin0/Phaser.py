#!/usr/bin/env python3

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import gzip
import random
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.SeqIO import FastaIO
except:
	print("Error: biopython module is not properly installed.")
	quit()

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

Phaser
Map reads, find SNPs, phase SNPs, extract consensus sequences

""")

########################################
############### ARGUMENTS ##############
########################################

parser.add_argument("-r","--reads",
					nargs='+',
					#default='sample1_R1.fastq.gz sample1_R2.fastq.gz',
					default=[],
					help="Space-separated list of reads in fastq(.gz) format (default: %(default)s)")
parser.add_argument("-t","--targets",
					type=str,
					default='../squamate_AHE_UCE_genes_loci2.fasta',
					help="Reference fasta to phase (default: %(default)s)")
parser.add_argument("-o","--output",
					type=str,
					default="09_phase",
					help="Folder in which to export results (default: %(default)s)")
parser.add_argument("-c","--cpu",
					type=int,
					default=8,
					help="Number of threads to be used in each step. (default: %(default)s)")
args=parser.parse_args()

########################################
################# SETUP ################
########################################

reads = [os.path.abspath(x) for x in args.reads]
targets_name = os.path.abspath(args.targets)

output = args.output
prefix = targets_name.split("/")[-1]
prefix = prefix.split(".")[:-1]
prefix = '.'.join(prefix)
prefix = os.path.join(output,prefix)
cpus = args.cpu

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting Phaser...")
start_dir = os.getcwd()
print("\tReads -> "+ ' '.join(args.reads))
print("\tTargets -> "+ args.targets)
print("\tOutput -> "+ output)
print("\tThreads -> "+ str(cpus))

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

def Unambigufy(seq):
	d = {'A': 'A',
		'C': 'C',
		'G': 'G',
		'T': 'T',
		'M': 'AC',
		'R': 'AG',
		'W': 'AT',
		'S': 'CG',
		'Y': 'CT',
		'K': 'GT',
		'V': 'ACG',
		'H': 'ACT',
		'D': 'AGT',
		'B': 'CGT',
		'-': 'N',
		'N': 'N'}
	
	seq2=''
	for letter in seq:
		choice=random.choice(d.get(letter))
		seq2= seq2 + choice
	return seq2

########################################
################# CODE #################
########################################

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Phasing ' + targets_name.split("/")[-1] + ' :::')

mkdir_p(output)

targets = list(SeqIO.parse(targets_name,'fasta'))

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Removing ambiguities prior to mapping :::')
# REMOVE AMBIGUITIES FOR PHASING
new_ref = []
for rec in targets:
	r = rec
	seq = str(r.seq).upper()
	new_seq = Unambigufy(seq)
	r.seq = Seq(new_seq)
	new_ref.append(r)

handle=open(prefix + '.tmp.fasta', "w")
writer = FastaIO.FastaWriter(handle)
writer.write_file(new_ref)
handle.close()

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Mapping and sorting reads :::')
sp.call('minimap2 -ax sr -t ' + str(cpus) + ' ' + prefix + '.tmp.fasta ' + ' '.join(reads) + ' | samtools sort -@ ' + str(cpus) + ' - > ' + prefix + '.tmp.bam', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
sp.call('samtools index ' + prefix + '.tmp.bam', shell=True)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Piling reads :::')
sp.call('bcftools mpileup -I -f ' + prefix + '.tmp.fasta ' + prefix + '.tmp.bam | bcftools call -mv -V indels -Oz > ' + prefix + '.tmp.vcf.gz', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
sp.call('bcftools index ' + prefix + '.tmp.vcf.gz', shell=True)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: WhatsHap Phasing :::')
sp.call('whatshap phase --ignore-read-groups --reference=' + prefix + '.tmp.fasta -o ' + prefix + '.tmp.phased.vcf ' + prefix + '.tmp.vcf.gz ' + prefix + '.tmp.bam', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
sp.call('bcftools norm -f ' + prefix + '.tmp.fasta -m +any -Oz -o ' + prefix + '.tmp.phased.norm.vcf.gz ' + prefix + '.tmp.phased.vcf', shell=True)
sp.call('tabix ' + prefix + '.tmp.phased.norm.vcf.gz', shell=True)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Extracting Consensus Fastas :::')
sp.call('bcftools consensus -H 1 -f ' + prefix + '.tmp.fasta ' + prefix + '.tmp.phased.norm.vcf.gz > ' + prefix + '.tmp.hap0.fasta', shell=True)
sp.call('bcftools consensus -H 2 -f ' + prefix + '.tmp.fasta ' + prefix + '.tmp.phased.norm.vcf.gz > ' + prefix + '.tmp.hap1.fasta', shell=True)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Combining Haplotypes :::')
combined_haps = []
hap0 = list(SeqIO.parse(prefix + '.tmp.hap0.fasta','fasta'))
hap1 = list(SeqIO.parse(prefix + '.tmp.hap1.fasta','fasta'))
for seq in hap0:
	seq.id = seq.name = seq.description = seq.id + "|h0"
	seq.seq = seq.seq.ungap("N")
	combined_haps.append(seq)

for seq in hap1:
	seq.id = seq.name = seq.description = seq.id + "|h1"
	seq.seq = seq.seq.ungap("N")
	combined_haps.append(seq)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Removing temporary files and writing final results :::')
sp.call('rm ' + prefix + '.tmp.*', shell=True)

handle=open(prefix + '.haplotypes.fasta', "w")
writer = FastaIO.FastaWriter(handle)
writer.write_file(combined_haps)
handle.close()
