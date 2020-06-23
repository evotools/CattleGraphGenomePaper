# Get unique sequences
## Introduction
The different script are run in sequence to generate a final bed coverage file, with the coverage support vector for each base in the genome, it extract the sequences and generate metrics on them.

## Dependencies
The pipeline requires the following software to be installed and accessible:
1. [hal](https://github.com/ComparativeGenomicsToolkit/hal)
2. [bedtools](https://bedtools.readthedocs.io/en/latest/)
3. [samtools](https://samtools.github.io/)
4. [CpGProD](http://doua.prabi.fr/software/cpgprod)
5. [RepeatMasker](http://www.repeatmasker.org/)
6. [HOMER](http://homer.ucsd.edu/homer/)
7. [ORFfinder](https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/)

## The pipeline
The pipeline consists of four main steps:
1. Conver the HAL archive to MAF by sequence to achieve speed and high parallelisation
2. Process the MAF to generate the different interval files.
3. Extract the sequences for each population-specific interval
4. Collect metrics on these sequences

## Running the pipeline
The scripts to run are numbered in the oreder to be run.
The scripts are meant to be used within a gridEngine environment, but can easily be 
adapted to other systems. 
To run the scripts, simply proceed as follow:
```
# Convert the hal file to maf for each sequence 
qsub -t 1-5 \
      -N hal2maf \
      00-HAL2MAF.sh hereford angus ndama ankole brahman BTAU.hal 

# Get coverage and metric for each genome 
qsub -t 1-5 \
      -N STATS \
      -hold_jid_ad hal2maf \
      01-GETINTERVALS.sh hereford angus ndama ankole brahman

# Process each genome, combining intervals if below 5 bp of distance and 
# keep if >100bp long 
qsub -t 1-5 \
      -N GETSEQ \
      -hold_jid_ad STATS \
      02-GETFASTA.sh hereford angus ndama ankole brahman 5 100 

# Characterise unique sequences in the genome
qsub -t 1-5 \
      -N SEQSTAT \
      -hold_jid_ad GETSEQ \
      03-GETSEQSTATS.sh hereford angus ndama ankole brahman
```

## Single script description
Below, we describe the output of each single script, which can easily run as standalone software
in case specific temporary files are needed (although some of can be *very large*).

### MAFtoCoverage
The script ```MAFtoCoverage.py``` gets maf files as inputs and generate a file with the following fields:
1. sequence
2. position
3. position
4. a hashtag-separed vector with:
   a) the number of alignments for each genome considered (e.g. ndama=0#hereford=2#angus=1)
   b) whether the base is masked and 
   c) whether is an N.
The script is parallelised and greatly benefits from it. 
The output is an unsorted bed file, although the file is not fully compliant at the moment.
To make it compliant, simply add 1 to the third position as follow:
```awk '{print $1,$2,$3+1,$4}' input.bed > output.bed```

### CoverageToSupport
```CoverageToSupport.py``` conver the data from the previous coverage file to a support vector format. 
This is a teb-separed file with the following fields:
1. sequence
2. position
3. position
4. support vector (e.g. if the genomes provided are hereford,angus,ndama, and there are alignments 
    between angus and ndama, the support vector will be 011).
5. whether is repetitive
6. whether is hard-masked (N)
As for the previous script, the file is a non-compliant bed file. To convert it to a compliant bed, just 
follow the instructions for ```MAFtoCoverage.py```.

### BigCoverageToSmallCoverage
```BigCoverageToSmallCoverage.py``` convert the large support to a concise (and compliant this time) bed file.
This is a simple bed file with the following fields:
1. sequence
2. initial position
3. final position
3. a hash-separed vector with:
   a) support vector (as described above)
   b) whether the base is masked (0/1 for unmasked/masked), 
   c) whether is hard-masked ((0/1 for not-N/N), and 
   d) number of bases in the region

It is possible and easy to extract unique, non repetitive regions by using combinations of support vectors, boolean filtering for 
repetitiveness and gaps and for size based on the windows size, in addition to filtering by position. 