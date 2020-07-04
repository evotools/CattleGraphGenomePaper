#!/bin/sh
. /etc/profile.d/modules.sh
module load anaconda
module load R
module load roslin/perl/5.26.1
module load roslin/samtools
module load roslin/bedtools
source activate DataPy

# Some pre-defined paths (less arguments from command line...)
export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/CpGProD/

# Get otions from command line.
while getopts ":f:o:b:n" opt; do
  case $opt in
    f) FASTA=${OPTARG};;
    o) OUTPUT=${OPTARG};;
    b) BEDFILE=${OPTARG};;
    n) NVERSION=${OPTARG};;
  esac
done

# Create temporary folder if not exists
if [ ! -e ${OUTPUT} ]; then
        mkdir ${OUTPUT}
fi

# Scanning for motifs into the genome
if [ ! -e ${OUTPUT}/CPG ]; then mkdir ${OUTPUT}/CPG; fi
CpGProD ${FASTA} \
        ${OUTPUT}/CPG/CPGPROD.out \
        -html ${OUTPUT}/CPG/CPGPROD.html
python -c "import sys;sys.stdout.write(''.join([line for line in open(sys.argv[1]) if len(line.strip().split()) > 2]) )" ${OUTPUT}/CPG/CPGPROD.out | \
  awk 'BEGIN{OFS="\t"};NR>1 {print $1,$3,$4,$5"#"$8}'> ${OUTPUT}/CPG/CPGPROD.bed
bedtools intersect -a ${OUTPUT}/CPG/CPGPROD.bed -b $BEDFILE -wa -wb > ${OUTPUT}/CPG/CPGPROD.specific.bed