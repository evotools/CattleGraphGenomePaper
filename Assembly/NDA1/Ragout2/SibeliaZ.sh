#!/bin/bash
# Run as follow:
# ./SibeliaZ.sh genome1.fa genome2.fa genome3.fa .. genomeN.fa 
# Each genome needs the name of the sequence to be in UCSC format (>genome.sequenceID)

# Combine all input genomes
cat ${@} 

# Run sibeliaZ
sibeliaz -f225 -t${NSLOTS} -o MyData sequences.fasta -k 21
