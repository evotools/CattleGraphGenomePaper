#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N graphLen
#$ -cwd
#$ -l h_rt=47:00:00
#$ -pe sharedmem 1
#$ -R y
#$ -l h_vmem=64.0G
#$ -hold_jid graphIds
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.err
#$ -P roslin_ctlgh
. /etc/profile.d/modules.sh
module load java/jdk/1.8.0
module load roslin/samtools/1.9
module load igmm/compilers/gcc

vg=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/vg


# 
# Start program
# 
export TMPDIR=`pwd`

for i in {1..29}; do $vg/vg stats -l GRAPH/CHR${i}.final.vg | awk -v var=${i} '{print var, $2}' ; done > pathLengths.txt

echo "Done paths length calculation"

