#!/bin/bash
#$ -cwd
#$ -N wtdbg2
#$ -pe sharedmem 8
#$ -l h_rt=96:00:00
#$ -l h_vmem=32G

# Prepare dependencies
. /etc/profile.d/modules.sh
module load roslin/samtools
module load roslin/bwa/0.7.17
wtdbg2F="/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/wtdbg2"
minimap="/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/minimap2"
ucscExe='/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/UCSC'
miniasm="/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/miniasm"

# Getting arguments
while getopts ":f:a:b:o:v" opt; do
  case $opt in
    f) inputFasta=${OPTARG};;
    a) read1=${OPTARG};;
    b) read2=${OPTARG};;
    o) outname=${OPTARG};;
    v) version=${OPTARG};;
  esac
done

# Running wtdbg2
${wtdbg2F}/wtdbg2 -x rs -g 2.7g -i $inputFasta -t ${NSLOTS} -fo $outname
# Derive consensus
${wtdbg2F}/wtpoa-cns -t ${NSLOTS} -i ${outname}.ctg.lay.gz -fo ${outname}.ctg.fa
echo "Done assembly, and start PacBio polishing."

