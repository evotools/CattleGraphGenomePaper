#!/bin/sh

. /etc/profile.d/modules.sh
module load anaconda
module load R
module load roslin/perl/5.26.1
module load roslin/samtools
module load roslin/bedtools
source activate DataPy

# Some pre-defined paths (less arguments from command line...)
export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/HOMER/bin/

# Get otions from command line.
while getopts ":b:g:o:n" opt; do
  case $opt in
    b) BED=${OPTARG};;
    g) GENOME=${OPTARG};;
    o) OUTPUT=${OPTARG};;
    n) NVERSION=${OPTARG};;
  esac
done

# Create temporary folder if not exists
if [ ! -e ./TMP ]; then
        mkdir ./TMP
fi
if [ ! -e ${OUTPUT}/MOTIFS ]; then
        mkdir ${OUTPUT}/MOTIFS
fi

# Scanning for motifs into the genome
findMotifsGenome.pl $BED $GENOME ${OUTPUT}/MOTIFS -mask -mset vertebrates 
scanMotifGenomeWide.pl /exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/HOMER/data/knownTFs/vertebrates/known.motifs \
                        $GENOME \
                        -bed -int -keepAll > ${OUTPUT}/MOTIFS/All.bed
bedtools intersect -a ${OUTPUT}/MOTIFS/All.bed -b $BED -wa >  ${OUTPUT}/MOTIFS/SPECIFIC.bed
