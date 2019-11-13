#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N graphJoin
#$ -cwd
#$ -l h_rt=8:00:00
#$ -pe sharedmem 1
#$ -R y
#$ -l h_vmem=8.0G
#$ -hold_jid graphAdd
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

cd GRAPH
cat $(for i in $(seq 1 29); do echo *CHR${i}.final.vg; done) > all.vg
echo "Done joint graph"

