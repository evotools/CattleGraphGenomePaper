#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N graphGCSA
#$ -cwd
#$ -l h_rt=160:00:00
#$ -pe sharedmem 12
#$ -R y
#$ -l h_vmem=80.0G
#$ -hold_jid graphAdd,graphIds,graphPrune,graphMap
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.err
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
${vg}/vg index -p -Z 4096 -t $NSLOTS -g all.gcsa -f mapping $(for i in $(seq 1 29); do echo chr${i}.final.pruned.vg; done) && echo "Done GCSA2"

