#!/bin/bash
#$ -cwd
#$ -N splitfa
#$ -pe sharedmem 1
#$ -R y
#$ -l h_rt=00:59:59
#$ -l h_vmem=16G
#$ -P roslin_ctlgh
#$ -o LOGS/splitfa.$TASK_ID.out
#$ -e LOGS/splitfa.$TASK_ID.err
. /etc/profile.d/modules.sh
module load igmm/apps/hmmer/3.1b2
module load roslin/samtools
module load roslin/blast+
module load igmm/apps/BEDOPS/2.4.26 

export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/RepeatMasker/Dependencies/bins
TRFPATH=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/


while getopts ":f:l:v" opt; do
  case $opt in
    f) fasta=${OPTARG};;
    l) length=${OPTARG};;
    v) version=${OPTARG};;
  esac
done

if [ ! -e ${fasta}.fai ]; then samtools faidx ${fasta}; fi
python SplitFastaEvenly.py ${fasta}.fai ${length}
mkdir SPLIT
for i in $(ls LIST/); do 
    fname=$(basename $i ".txt" )
    samtools faidx -r LIST/$i ndama.fasta > SPLIT/${fname}.fasta
    echo SPLIT/${fname}.fasta
done > fastaList.txt