#!/bin/bash
# Script used to generate the linear expanded graph genome.
# Grid Engine options (lines prefixed with #$ or #!)
#$ -N graphAdd
#$ -cwd
#$ -l h_rt=47:00:00
#$ -pe sharedmem 1
#$ -R y
#$ -l h_vmem=48.0G
#$ -hold_jid gc_60419
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
. /etc/profile.d/modules.sh
module load java/jdk/1.8.0
module load roslin/samtools/1.9
module load igmm/compilers/gcc
module load igmm/apps/tabix
export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/vg


# 
# Start program
# 
while getopts ":f:v:o:n" opt; do
  case $opt in
    f) fa=${OPTARG};;
    v) vcf=${OPTARG};;
	o) oname=${OPTARG};;
    n) nversion=${OPTARG};;
  esac
done




export TMPDIR=`pwd`

echo $vcf
echo $fa

vg construct -p -v $vcf -r ${fa} > $oname.vg