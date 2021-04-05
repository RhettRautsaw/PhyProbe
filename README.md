# PhyProbe

## Pipeline

[Flowchart](https://docs.google.com/drawings/d/1YbIQsHFdcLFzGP--4vbma8IAR7dfUapjAgoFcXb25_E/edit?usp=sharing)

1. Assembly: rnaSpades or Trinity
2. BUSCO 4.1.4: tetrapoda_odb10
	- BuscoCleaner.py:
		- Check for duplicates in multi-copy and move to single copy
		- Concatenate all individuals into busco_loci.fasta
3. CodAn: VERT_full
4. Extract longest isoforms (LongForms.py)
5. cd-hit-est: 100% 
6. Extract sqCL loci (LocusExtractor.py)
7. Orthofinder on clustered ORFs
	- scOrthoSorter.py:
		- Extract single copy loci with less than 50% missing individuals
8. Combine Data Sources
	- UniqRenamer.py:
		- Compare BUSCO loci to sqCL targets.
		- Extract matching regions and rename appropriately
		- This makes sure that we don't have repeat loci for phylogenetics later
	- UniqRenamer.py:
		- Compare Orthofinder loci to combined sqCL and BUSCO loci
	- SampleSorter.py:
		- Sort combined BUSCO, OF, and sqCL extracted loci into each individual
9. Phasing (Phaser.py)

*** 

## Installation
```
python3
	- biopython
	- pyfastx
	- pandas
	- numpy
	- dfply
	- p_tdqm
	- tdqm
	- multiprocessing
	- glob

rnaSpades
BUSCO
CodAn
BLAST
Minimap2
cd-hit
seqtk
samtools
bcftools
WhatsHap
orthofinder

# Conda Install
perl
perl-bioperl
perl-mce


```

## Directory Setup
```
PhyProbe_Directory/
├── Other
|	├── sample_001							paired-end read sample
|	|	├── sample_001_R1.fastq.gz
|	|	└── sample_001_R2.fastq.gz
|	└── sample_002							single-end read sample
|		└── sample_002_R1.fastq.gz
├── SeqCap
|	├── sample_001
|	|	├── sample_001_R1.fastq.gz
|	|	└── sample_001_R2.fastq.gz
|	└── sample_002
|		└── sample_002_R1.fastq.gz
├── targets.fasta
└── Transcriptomes
	├── sample_001
	|	├── sample_001_R1.fastq.gz
	|	└── sample_001_R2.fastq.gz
	└── sample_002
		└── sample_002_R1.fastq.gz
```