#!/usr/bin/env python

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import multiprocessing as mp
from psutil import virtual_memory
import glob

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

       __   ____  _           ____            _           
   ___|    |  _ \| |__  _   _|  _ \ _ __ ___ | |__   ___  
  |   |__  | |_) | '_ \| | | | |_) | '__/ _ \| '_ \ / _ \ 
--|        |  __/| | | | |_| |  __/| | | (_) | |_) |  __/ 
  |______  |_|   |_| |_|\__, |_|   |_|  \___/|_.__/ \___| 
                        |___/                             

PhyProbe is a bioinformatic pipeline designed to extract phylogenetic loci from Next-Generation 
Sequencing datasets including RNAseq, WGS, and Sequence/Target Capture methodologies (e.g., AHE, UCEs).
Specifically, PhyProbe extracts three categories of loci:

	1. BUSCO Loci
	2. Transcripts (RNAseq data only)
	3. Provided Target Loci (e.g., AHE or UCEs)

If using multiple datatypes (multi-omics), we recommend processing RNAseq data first since it is a very
versatile datatype. By processing RNAseq data first, you can capture BUSCOs and Transcripts that can be used as 
targets for the other datatypes. Not only will this speed up processing, but ensures that paralogs are removed
from BUSCO results and that additional Transcript-based loci can be captured from other datatypes.

#### PIPELINE ####

[ Individual-sample analysis modes ]
-m rna: assembly (rnaspades), BUSCO, transcript capture (CodAn/OrthoFinder), target capture (AliBaSeq), clustering (cd-hit-est), phasing (WhatsHap)
-m wgs: assembly (MaSuRCA), BUSCO, target capture, clustering, phasing
-m other: assembly (rnaspades), BUSCO, target capture, clustering, phasing

[ Grouped-sample analysis modes ]
-m ortho: OrthoFinder
-m mcl: MCL
-m phylo: combine samples, separate loci, align (mafft), trim (trimal & CIAlign), outlier removal (IQTree & Treeshrink), missing data filtering, & gene tree inference (IQTree)
""")

parser.add_argument("-m", "--mode",
					type=str,
					default=None,
					help="Analysis mode. Single sample options: [rna, wgs, other]. All sample options: [mcl, ortho, phylo] (default: %(default)s)")
parser.add_argument("-i","--input",
					type=str,
					default=None,
					help="For -m [rna, wgs, other]: a folder with DNA/RNA read data and/or previously completed steps for a single sample. Folder name should correspond to sample name as it will used as a prefix in each step. \n For -m [mcl, ortho, phylo]: a txt file with a list of samples/folders. PhyProbe will autodetect which step to resume a previous analysis. (default: %(default)s)")
parser.add_argument("-d","--delim",
					type=str,
					default='|',
					help="Delimiter used to separate gene or locus names from sample names. For example, Locus0001|SampleName. (default: %(default)s)")

group1=parser.add_argument_group('BUSCO', 'Options for BUSCO locus capture.')
group1.add_argument("--busco",
					type=str,
					default=None,
					help="Required for BUSCO step. Path to downloaded BUSCO db, name of BUSCO db to download, or \'auto\'. If not provided identification of BUSCO loci will be skipped (unless already finished and provided in target database). (default: %(default)s)")

group3=parser.add_argument_group('Transcript Capture', 'Options for transcript locus capture from RNAseq.')
group3.add_argument("--codan",
					type=str,
					default=None,
					help="Required for CodAn step. Path to CodAn model for transcriptome CDS prediction. If not provided identification of additional transcripts will be skipped (unless already finished and provided in target database). (default: %(default)s)")

group2=parser.add_argument_group('Target Capture', 'Options for target locus capture (from provided reference) with AliBaSeq.')
group2.add_argument("--targets",
					type=str,
					default=None,
					help="Required for AliBaSeq step. Fasta with targets from which to extract loci. If not provided identification of target loci will be skipped. (default: %(default)s)")

group4=parser.add_argument_group('Processing', 'Options for speeding up processing steps.')
group4.add_argument("-c","--cpu",
					type=int,
					default=mp.cpu_count(),
					help="Number of threads to be used in each step. (default: %(default)s)")
group4.add_argument("-mem","--memory",
					type=int,
					default=(round(virtual_memory().total/1e9)-1),
					help="Max amount of memory for assembly. Only applicable for rnaSpades assembly. (default: %(default)s)")
group4.add_argument("--nodes",
					type=str,
					default=None,
					help="File with list of HPC nodes that you can ssh into to parallelize alignments and trimming. Only applicable for -m phylo. (default: %(default)s)")
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

def rnaspades(reads, cpus, mem):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: rnaSpades Assembly :::')
	if len(reads)==1:
		sp.call("rnaspades.py -s " + reads[0] + " -t " + str(cpus) + " -m " + str(mem) + " -o 01_assembly", shell=True, stdout=errlog, stderr=errlog)
	elif len(reads)==2:
		sp.call("rnaspades.py -1 " + reads[0] + " -2 " + reads[1] + " -t " + str(cpus) + " -m " + str(mem) + " -o 01_assembly", shell=True, stdout=errlog, stderr=errlog)
	elif len(reads)>2 or len(reads)==0:
		print(str(len(reads)) + " fastq files identified")
		print("Please only include one file for single-end mode or two files for paired-end mode")
		quit()
	sp.call("mv 01_assembly/transcripts.fasta 01_assembly.fasta", shell=True, stdout=errlog, stderr=errlog)

def masurca(reads, cpus):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: MaSuRCA :::')
	sp.call("mkdir 01_assembly", shell=True, stdout=errlog, stderr=errlog)
	os.chdir("01_assembly")
	sp.call("masurca -t " + str(cpus) + " -i " + ','.join(reads), shell=True, stdout=errlog, stderr=errlog)
	sp.call("mv CA/primary.genome.scf.fasta " + input + "/01_assembly.fasta", shell=True, stdout=errlog, stderr=errlog)
	os.chdir(input)

def busco(assembly, busco_mod, cpus, prefix):
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: BUSCO-' + busco_mod + ' :::')
	sp.call("busco -i " + assembly + " -c " + str(cpus) + " -o 02_busco -m geno -l " + busco_mod, shell=True, stdout=errlog, stderr=errlog)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Cleaning BUSCO Results :::')
	sp.call("BuscoCleaner.py -f 02_busco -a " + assembly + " -n " + prefix + " -d \'" + delim + "\' -c " + str(cpus), shell=True, stdout=errlog, stderr=errlog)

def codan(assembly, codan_mod, cpus, prefix, delim):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Protein-Coding Gene Identification :::')
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: CodAn :::')
	sp.call("codan.py -t " + assembly + " -m " + codan_mod + " -o  03_codan -c " + str(cpus), shell=True, stdout=errlog, stderr=errlog)
	os.chdir("03_codan")
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Extracting Longest Isoforms :::')
	sp.call("LongForms.py -i ORF_sequences.fasta -o ORF_sequences.long.fasta", shell=True, stdout=errlog, stderr=errlog)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Clustering @ 100% Identity :::')
	sp.call("cd-hit-est -i ORF_sequences.long.fasta -o ORF_sequences.clust.fasta -c 1 -d 0", shell=True, stdout=errlog, stderr=errlog)
	sp.call("perl -pi -e 's/>/>" + prefix + delim +"/g' ORF_sequences.clust.fasta", shell=True, stdout=errlog, stderr=errlog)
	os.chdir(input)

def alibaseq(assembly, targets, cpus, prefix, delim, output):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Target Capture :::')
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Make BLAST Database :::')
	sp.call("mkdir " + output + "; cp " + assembly + " " + output + "/tmp_assembly.fasta", shell=True, stdout=errlog, stderr=errlog)
	os.chdir(output)
	sp.call("makeblastdb -in tmp_assembly.fasta -dbtype nucl -parse_seqids", shell=True, stdout=errlog, stderr=errlog)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: BLAST :::')
	sp.call("blastn -task dc-megablast -query " + targets + " -db tmp_assembly.fasta -outfmt 6 -out tmp_assembly.fasta.blast -evalue 1e-10 -num_threads " + str(cpus), shell=True, stdout=errlog, stderr=errlog)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: AliBaSeq :::')
	sp.call("alibaseqPy3.py -x a -f S -b tmp_assembly.fasta.blast -t tmp_assembly.fasta -e 1e-10 --is --amalgamate-hits", shell=True, stdout=errlog, stderr=errlog)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Cleaning AliBaSeq Results :::')
	sp.call("AliBaSeqCleaner.py -f alibaseq_out -r " + targets + " -n " + prefix + " -d \'" + delim + "\' -c " + str(cpus), shell=True, stdout=errlog, stderr=errlog)
	sp.call("rm -r alibaseq_out *.log tmp_*", shell=True, stdout=errlog, stderr=errlog)
	sp.call("mv alibaseq_clean/* .", shell=True, stdout=errlog, stderr=errlog)
	sp.call("rm -r alibaseq_clean", shell=True, stdout=errlog, stderr=errlog)
	os.chdir(input)

########################################
################# SETUP ################
########################################

if args.input==None:
	print("No input specified. Use -i/--input flag to specify input")
	quit()
else:
	input = os.path.abspath(args.input)
	prefix = os.path.basename(input)

if args.mode==None:
	print("No mode specified. Use -m/--mode flag to specify \'DNA\', \'RNA\', or \'GROUP\' depending on which step of the PhyProbe pipeline you are at")
	quit()
else:
	mode = args.mode
	if mode not in ['rna', 'wgs', 'other', 'mcl', 'ortho', 'phylo']:
		print("Mode specified did not match any of the possible options [rna, wgs, other, mcl, ortho, phylo]")
		quit()
	if mode in ['rna', 'wgs', 'other'] and not os.path.isdir(input):
		print("Mode was specified as -m " + mode + ", but -i is not a directory.")
		quit()
	if mode in ['mcl', 'ortho', 'phylo'] and not os.path.isfile(input):
		print("Mode was specified as -m " + mode + ", but -i is not a file")
		quit()

delim = args.delim
cpus = args.cpu
mem = args.memory

errlog = open("PhyProbe.errlog", 'a')

########################################
############### PREAMBLE ###############
########################################

print("""
       __   ____  _           ____            _           
   ___|    |  _ \| |__  _   _|  _ \ _ __ ___ | |__   ___  
  |   |__  | |_) | '_ \| | | | |_) | '__/ _ \| '_ \ / _ \ 
--|        |  __/| | | | |_| |  __/| | | (_) | |_) |  __/ 
  |______  |_|   |_| |_|\__, |_|   |_|  \___/|_.__/ \___| 
                        |___/                             
""")

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting PhyProbe v1.0 :::")
print("\tMode:\t"+ mode)
print("\tInput:\t"+ input)
print("\tDelim:\t" + delim)
if args.busco==None:
	busco_mod=None
	print("\tBUSCO:\tSkipping")
else:
	busco_mod = os.path.abspath(args.busco)
	print("\tBUSCO:\t" + busco_mod)
if args.codan==None:
	codan_mod=None
	print("\tCodAn:\tSkipping")
else:
	codan_mod = os.path.abspath(args.codan)
	print("\tCodAn:\t" + codan_mod)
if args.targets==None:
	targets=None
	print("\tAliBaSeq:\tSkipping")
else:
	targets = os.path.abspath(args.targets)
	print("\tAliBaSeq:\t" + targets)
print("\tCPUs:\t" + str(cpus))
print("\tMem:\t" + str(mem))

########################################
################# START ################
########################################
start_dir = os.getcwd()

if mode in ['rna', 'wgs', 'other']:
	os.chdir(input)
	reads = [os.path.abspath(x) for x in glob.glob("*.fastq*")+ glob.glob("*.fq*")]
	reads.sort()
	
	# O1_assembly
	if not os.path.isfile("01_assembly.fasta"):
		if mode in ['rna', 'other']:
			rnaspades(reads, cpus, mem)
		else:
			masurca(reads, cpus)
	
	sp.call("cp 01_assembly.fasta 01_renamed.fasta", shell=True, stdout=errlog, stderr=errlog)
	sp.call("FastaRenamer.py -f 01_renamed.fasta -n contig", shell=True, stdout=errlog, stderr=errlog)
	
	# 02_busco
	if not os.path.isfile("02_busco/" + prefix + "_busco.fasta") and busco_mod != None:
		if os.path.isfile(start_dir + "/00_mcl/BUSCO_db.fasta"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: MCL results detected! :::')
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Using AliBaSeq to extract orthologous BUSCOs :::')
			alibaseq("01_renamed.fasta", start_dir + "/00_mcl/BUSCO_db.fasta", cpus, prefix, delim, "02_busco")
			sp.call("mv " + prefix + "_targets.fasta " + prefix + "_busco.fasta", shell=True, stdout=errlog, stderr=errlog)
		else:
			busco("01_renamed.fasta", busco_mod, cpus, prefix, delim)
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: WARNING: We\'ve noticed that BUSCO extraction results in a lot of sequences that are likely paralogs and will impact downstream alignment and phylogenetic analyses. Once all your samples finish with BUSCO, we highly recommend that you run \'-m mcl\' with all your samples to cluster sequences. We can then keep the largest cluster as our reference ortholog and re-extract more confident orthologs from the rest with AliBaSeq. :::')
	
	# 03_codan
	if not os.path.isfile("03_codan/" + prefix + "_transcripts.fasta") and codan_mod != None:
		if os.path.isfile(start_dir + "/00_ortho_res/OrthoCleaner_renamed.fasta"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: OrthoFinder results detected! :::')
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Using AliBaSeq to extract single-copy transcript Orthogroups :::')
			alibaseq("01_renamed.fasta", start_dir + "/00_ortho_res/OrthoCleaner_renamed.fasta", cpus, prefix, delim, "03_codan")
			sp.call("mv " + prefix + "_targets.fasta " + prefix + "_transcripts.fasta", shell=True, stdout=errlog, stderr=errlog)
		elif not os.path.isfile("03_codan/ORF_sequences.clust.fasta") and mode == 'rna':
			codan("01_assembly.fasta", codan_mod, cpus, prefix, delim)
		if not os.path.isfile(start_dir + "/00_ortho_res/OrthoCleaner_renamed.fasta"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: No OrthoFinder results detected..yet. I would guess that your next step is to detect single-copy orthologs across all your samples. To do this, once all your RNA samples complete the CodAn step, run \'-m ortho\' with a list of all your samples/directories. Once complete, come back and rerun this same command to re-extract protein-coding loci with AliBaSeq. :::')
			print('\n' + dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") +' ::: PS...Now might be a good time to run \'-m mcl\' mode too if you are looking for BUSCO loci! You can then add the clustered BUSCO loci to the target database and re-extract! :::')
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Quitting PhyProbe now... :::')
			quit()
	
	# 04_targets
	if not os.path.isfile("04_targets/" + prefix + "_targets.fasta") and targets != None:
		alibaseq("01_renamed.fasta", targets, cpus, prefix, delim, "04_targets")
	
	# 05_combining
	if not os.path.isfile("05_combine/" + prefix + "_nonredun.fasta"):
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Removing Redundancy Across BUSCO, Target, and Transcript Results :::')
		mkdir_p("05_combine")
		if os.path.isfile("02_busco/" + prefix + "_targets.fasta"):
			sp.call("Ungap.py -f 02_busco/" + prefix + "_targets.fasta -o 05_combine/" + prefix + "_busco.fasta", shell=True, stdout=errlog, stderr=errlog)
		else:
			sp.call("Ungap.py -f 02_busco/" + prefix + "_busco.fasta -o 05_combine/" + prefix + "_busco.fasta", shell=True, stdout=errlog, stderr=errlog)
		sp.call("Ungap.py -f 03_codan/" + prefix + "_targets.fasta -o 05_combine/" + prefix + "_transcripts.fasta", shell=True, stdout=errlog, stderr=errlog)
		sp.call("Ungap.py -f 04_targets/" + prefix + "_targets.fasta -o 05_combine/" + prefix + "_targets.fasta", shell=True, stdout=errlog, stderr=errlog)
		os.chdir("05_combine")
		# Find Targets not already represented by BUSCOs and concatenate
		sp.call("cd-hit-est-2d -i " + prefix + "_busco.fasta -i2 " + prefix + "_targets.fasta -o " + prefix + "_targets2.fasta -c 0.99 -d 0", shell=True, stdout=errlog, stderr=errlog)
		sp.call("cat " + prefix + "_busco.fasta " + prefix + "_targets2.fasta > " + prefix + "_combo1.fasta", shell=True, stdout=errlog, stderr=errlog)
		# Find Transcripts not already represented by BUSCOs or Targets and concatenate
		sp.call("cd-hit-est-2d -i " + prefix + "_combo1.fasta -i2 " + prefix + "_transcripts.fasta -o " + prefix + "_transcripts2.fasta -c 0.99 -d 0", shell=True, stdout=errlog, stderr=errlog)
		sp.call("cat " + prefix + "_combo1.fasta " + prefix + "_transcripts2.fasta  > " + prefix + "_nonredun.fasta", shell=True, stdout=errlog, stderr=errlog)
		os.chdir(input)
	
	# 06_phasing
	if not os.path.isfile("06_phased/" + prefix + "_phased.fasta"):
		mkdir_p("06_phased")
		sp.call("cp 05_combine/" + prefix + "_nonredun.fasta 06_phased/", shell=True, stdout=errlog, stderr=errlog)
		os.chdir("06_phased")
		sp.call("VariantCaller.py -f " + ' '.join(reads) + " -r " + prefix + "_nonredun.fasta -s " + prefix + " -c " + str(cpus) + " --mpileup --haplotypes", shell=True, stdout=errlog, stderr=errlog)
		sp.call("mv 11_haplotmp1.fasta 11_haplotype1.fasta; mv 11_haplotmp2.fasta 11_haplotype2.fasta", shell=True, stdout=errlog, stderr=errlog)
		sp.call("cat 11_haplotype1.fasta 11_haplotype2.fasta > " + prefix + "_phased.fasta", shell=True, stdout=errlog, stderr=errlog)
		sp.call("Ungap.py -f " + prefix + "_phased.fasta", shell=True, stdout=errlog, stderr=errlog)

############################
###### GROUP ANALYSES ######
############################

if mode in ['mcl', 'ortho', 'phylo']:
	# Read List
	samples = open(input, "r").read().split("\n")
	samples = [i for i in samples if i]
	
	# 00_ortho
	if mode == 'ortho':
		mkdir_p("00_ortho")
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Copying CodAn results for all samples :::')
		for prefix in samples:
			if os.path.isfile(prefix + "/03_codan/ORF_sequences.clust.fasta"):
				sp.call("cp " + prefix + "/03_codan/ORF_sequences.clust.fasta 00_ortho/" + prefix + "_ORFs.fasta", shell=True, stdout=errlog, stderr=errlog)
			else:
				print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: No CodAn results detected for ' +  prefix + ' :::')
		# Running OrthoFinder
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Running OrthoFinder :::')
		sp.call("orthofinder -f 00_ortho -t " + str(cpus) + " -d -og -o 00_ortho_res", shell=True, stdout=errlog, stderr=errlog)
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Extracting single-copy Orthogroups with less than 50% missing data :::')
		sp.call("OrthoCleaner.py -f 00_ortho_res -d \'" + delim + "\' -c " + str(cpus) + " -p 50", shell=True, stdout=errlog, stderr=errlog)
		print('\n' + dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") +' ::: FINISHED \'-m ortho\' mode. You can now rerun \'-m rna\' with \'--codan\' specified and PhyProbe will auto-detect the OrthoCleaner results and output the final results in the codan directory as *_targets.fasta ... OR ... you can concatenate the OrthoCleaner_renamed.fasta to your target database for all your samples and only supply \'--targets\'. Although the latter method will make combining your sequences a little more challenging later! Good luck! :::')
		quit()
	
	# 00_mcl
	if mode == 'mcl':
		mkdir_p("00_mcl")
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Concatenating BUSCO results for all samples :::')
		for prefix in samples:
			if os.path.isfile(prefix + "/02_busco/" + prefix + "_busco.fasta"):
				sp.call("cat " + prefix + "/02_busco/" + prefix + "_busco.fasta >> 00_mcl/all_sample_buscos.fasta", shell=True, stdout=errlog, stderr=errlog)
			else:
				print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: No BUSCO results detected for ' +  prefix + ' :::')
		# Separating Loci
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Separating Loci :::')
		os.chdir("00_mcl")
		sp.call("SeqSorter.py -i all_sample_buscos.fasta -d \'" + delim + "\' -p 0 -c " + str(cpus), shell=True, stdout=errlog, stderr=errlog)
		sp.call("ls SeqSorter/*.fasta | perl -pe 's/.*\///g' > 00_loci.list", shell=True, stdout=errlog, stderr=errlog)
		# MCL Cluster each Locus
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: MCL Clustering & Keeping Largest Cluster as Reference :::')
		sp.call("parallel -a 00_loci.list -j " + str(cpus) + " --bar ' MCL.py -i SeqSorter/{} --inflation 3 -m 1 -c 1 &> /dev/null'", shell=True, stdout=errlog, stderr=errlog)
		# New BUSCO Database
		sp.call("cat MCL_out/*.fasta > BUSCO_db.fasta", shell=True, stdout=errlog, stderr=errlog)
		sp.call("perl -pe 's/\\" + delim + ".*//g' BUSCO_db.fasta > BUSCO_db_renamed.fasta", shell=True, stdout=errlog, stderr=errlog)
		print('\n' + dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") +' ::: FINISHED \'-m mcl\' mode. You can now rerun \'-m [rna, wgs, other]\' with \'--busco\' specified and PhyProbe will auto-detect the MCL results and output the final results in the busco directory as *_targets.fasta ... OR ... you can concatenate the BUSCO_db_renamed.fasta to your target database for all your samples and only supply \'--targets\'. Although the latter method will make combining your sequences a little more challenging later! Good luck! :::')
		quit()
	
	# 01_phylo
	if mode == 'phylo':
		mkdir_p('01_phylo')
		if not os.path.isfile("01_phylo/all_sample_loci.fasta"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Concatenating final sample results :::')
			for prefix in samples:
				if os.path.isfile(prefix + "/06_phased/" + prefix + "_phased.fasta"):
					sp.call("cat " + prefix + "/06_phased/" + prefix + "_phased.fasta >> 01_phylo/all_sample_loci.fasta", shell=True, stdout=errlog, stderr=errlog)
				else:
					print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: No Phased Haplotypes detected for ' +  prefix + ' :::')
		
		os.chdir("01_phylo")
		if not os.path.isdir("SeqSorter") or not os.path.isfile("00_loci.list"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Sorting Sequences by Locus :::')
			sp.call("SeqSorter.py -i all_sample_loci.fasta -d \'" + delim + "\' -p 0 -c " + str(cpus), shell=True, stdout=errlog, stderr=errlog)
			sp.call("ls SeqSorter/*.fasta | perl -pe 's/SeqSorter\///g' | perl -pe 's/.fasta//g' > 00_loci.list", shell=True, stdout=errlog, stderr=errlog)
		
		#sp.call("mkdir 00_seqs 01_aln 02_CIAlign 03_trimal 04_dropgappy 05_aln 06_treeshrink", shell=True, stdout=errlog, stderr=errlog)
		
		# 00_seqs
		if not os.path.isdir("00_seqs"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Removing Locus Information from Fasta Headers :::')
			mkdir_p("00_seqs")
			if args.nodes != None:
				sp.call("parallel -a 00_loci.list -j 1 --workdir $PWD --bar --sshloginfile " + args.nodes + " \"perl -pe 's/^>.*?" + delim + "/>/g' SeqSorter/{}.fasta > 00_seqs/{}.fasta\"", shell=True, stdout=errlog, stderr=errlog)
			else:
				sp.call("parallel -a 00_loci.list -j " + str(cpus) + " --bar \"perl -pe 's/^>.*?\\" + delim + "/>/g' SeqSorter/{}.fasta > 00_seqs/{}.fasta\"", shell=True, stdout=errlog, stderr=errlog)
		
		# 01_aln
		if not os.path.isdir("01_aln"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: MAFFT Alignment :::')
			mkdir_p("01_aln")
			if args.nodes != None:
				sp.call("parallel -a 00_loci.list -j 1 --workdir $PWD --bar --sshloginfile " + args.nodes + " \"mafft --auto --adjustdirection accurately --thread " + str(cpus) + " 00_seqs/{}.fasta > 01_aln/{}.fasta; perl -pi -e 's/_R_//g' 01_aln/{}.fasta\"", shell=True, stdout=errlog, stderr=errlog)
			else:
				sp.call("parallel -a 00_loci.list -j " + str(cpus) + " --bar \"mafft --auto --adjustdirectionaccurately --thread 1 00_seqs/{}.fasta > 01_aln/{}.fasta; perl -pi -e 's/_R_//g' 01_aln/{}.fasta\"", shell=True, stdout=errlog, stderr=errlog)
		
		# 02_CIAlign
		if not os.path.isdir("02_CIAlign"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: CIAlign Cleaning :::')
			mkdir_p("02_CIAlign")
			if args.nodes != None:
				sp.call("parallel -a 00_loci.list -j 1 --workdir $PWD --bar --sshloginfile " + args.nodes + " \"CIAlign --infile 01_aln/{}.fasta --outfile_stem 02_CIAlign/{} --remove_divergent --remove_divergent_minperc 0.80 --remove_insertions --crop_ends --remove_short\"", shell=True, stdout=errlog, stderr=errlog)
			else:
				sp.call("parallel -a 00_loci.list -j " + str(cpus) + " --bar \"CIAlign --infile 01_aln/{}.fasta --outfile_stem 02_CIAlign/{} --remove_divergent --remove_divergent_minperc 0.80 --remove_insertions --crop_ends --remove_short\"", shell=True, stdout=errlog, stderr=errlog)
		
		# 03_trimal
		if not os.path.isdir("03_trimal"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: trimal gappyout Cleaning :::')
			mkdir_p("03_trimal")
			if args.nodes != None:
				sp.call("parallel -a 00_loci.list -j 1 --workdir $PWD --bar --sshloginfile " + args.nodes + "\"trimal -in 02_CIAlign/{}_cleaned.fasta -out 03_trimal/{}.fasta -gappyout\"", shell=True, stdout=errlog, stderr=errlog)
			else:
				sp.call("parallel -a 00_loci.list -j " + str(cpus) + " --bar \"trimal -in 02_CIAlign/{}_cleaned.fasta -out 03_trimal/{}.fasta -gappyout\"", shell=True, stdout=errlog, stderr=errlog)
		
		# 04_dropgappy
		if not os.path.isdir("04_dropgappy"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Dropping Gappy Samples from Alignments :::')
			sp.call("ls 03_trimal/*.fasta | perl -pe 's/03_trimal\///g' | perl -pe 's/.fasta//g' > 03_loci.list", shell=True, stdout=errlog, stderr=errlog)
			mkdir_p("04_dropgappy")
			if args.nodes != None:
				sp.call("parallel -a 03_loci.list -j 1 --workdir $PWD --bar --sshloginfile " + args.nodes + " \"DropGappy.py -i 03_trimal/{}.fasta -o 04_dropgappy/{}.fasta -c 70 -ug\"", shell=True, stdout=errlog, stderr=errlog)
			else:
				sp.call("parallel -a 03_loci.list -j " + str(cpus) + " --bar \"DropGappy.py -i 03_trimal/{}.fasta -o 04_dropgappy/{}.fasta -c 70 -ug\"", shell=True, stdout=errlog, stderr=errlog)
		
		# 05_trees
		if not os.path.isdir("05_trees"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Realigning (MAFFT) and Building Preliminary Trees (IQTree) :::')
			mkdir_p("05_trees")
			mkdir_p("06_treeshrink")
			if args.nodes != None:
				sp.call("parallel -a 03_loci.list -j 1 --workdir $PWD --bar --sshloginfile " + args.nodes + " " +
				"\"mkdir 06_treeshrink/{}; " + 
				"mafft --auto --thread " + str(cpus) + " 04_dropgappy/{}.fasta > 05_trees/{}.fasta; " + 
				"cd 05_trees; " + 
				"iqtree -s {}.fasta -B 1000 -T " + str(cpus) + "; " + 
				"cp {}.fasta ../06_treeshrink/{}/input.fasta; " + 
				"cp {}.fasta.treefile ../06_treeshrink/{}/input.tree\"", 
				shell=True, stdout=errlog, stderr=errlog)
			else:
				sp.call("parallel -a 03_loci.list -j " + str(cpus) + " --bar " + 
				"\"mkdir 06_treeshrink/{}; " + 
				"mafft --auto 04_dropgappy/{}.fasta > 05_trees/{}.fasta; " + 
				"cd 05_trees; " + 
				"iqtree -s {}.fasta -B 1000; " + 
				"cp {}.fasta ../06_treeshrink/{}/input.fasta; " +
				"cp {}.fasta.treefile ../06_treeshrink/{}/input.tree\"",
				shell=True, stdout=errlog, stderr=errlog)
				sp.call("touch 06_treeshrink/COMPLETE", shell=True, stdout=errlog, stderr=errlog)
		
		# 06_treeshrink
		if not os.path.isfile("06_treeshrink/COMPLETE"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Checking/Removing Branch Outliers (TreeShrink) :::')
			sp.call("run_treeshrink.py -i 06_treeshrink -t input.tree -a input.fasta -c", shell=True, stdout=errlog, stderr=errlog)
		
		# 07_trimal
		if not os.path.isdir("07_trimal"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: trimal gappyout Cleaning :::')
			mkdir_p("07_trimal")
			if args.nodes != None:
				sp.call("parallel -a 00_loci.list -j 1 --workdir $PWD --bar --sshloginfile " + args.nodes + " \"trimal -in 06_treeshrink/{}/output.fasta -out 07_trimal/{}.fasta -gappyout\"", shell=True, stdout=errlog, stderr=errlog)
			else:
				sp.call("parallel -a 00_loci.list -j " + str(cpus) + " --bar \"trimal -in 06_treeshrink/{}/output.fasta -out 07_trimal/{}.fasta -gappyout\"", shell=True, stdout=errlog, stderr=errlog)
		
		# 08_MDF
		if not os.path.isdir("08_subsets"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Filtering Missing Data :::')
			mkdir_p("08_subsets")
			sp.call("FilterMissData.py -f 07_trimal -o 08_subsets -c " + str(cpus), shell=True, stdout=errlog, stderr=errlog)
		
		# 09_GeneTrees
		if not os.path.isfile("09_genetrees.tre"):
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Using genes with > 25% total taxa and taxa with > 5% total genes as final dataset :::')
			print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ::: Inferring Final Gene Trees :::')
			sp.call("ls 08_subsets/Genes25Taxa_Taxa5Genes/*.fasta | perl -pe \'s/.*\///g' > 09_loci.list", shell=True, stdout=errlog, stderr=errlog)
			if args.nodes != None:
				sp.call("parallel -a 09_loci.list -j 1 --workdir $PWD --bar --sshloginfile " + args.nodes + " \"cd 08_subsets/Genes25Taxa_Taxa5Genes; " + 
				"iqtree -s {} -B 1000 -T " + str(cpus) + "\"", shell=True, stdout=errlog, stderr=errlog)
			else:
				sp.call("parallel -a 09_loci.list -j " + str(cpus) + " --bar \"cd 08_subsets/Genes25Taxa_Taxa5Genes; " + 
				"iqtree -s {} -B 1000 -T 1\"", shell=True, stdout=errlog, stderr=errlog)
			sp.call("cat 08_subsets/Genes25Taxa_Taxa5Genes/*.treefile > 09_genetrees.tre", shell=True, stdout=errlog, stderr=errlog)
