#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N graphViz
#$ -cwd
#$ -l h_rt=8:00:00
#$ -pe sharedmem 1
#$ -R y
#$ -l h_vmem=128.0G
#$ -hold_jid gc_60419
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.err
#$ -P roslin_ctlgh
. /etc/profile.d/modules.sh
module load java/jdk/1.8.0
module load roslin/samtools/1.9

vg=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/vg
gviz=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/GVIZ/bin/

# 
# Start program
# 

bname=`basename -s ".xg" $1`

${vg}/vg viz -x $1 -o ${bname}.full.dot
${vg}/vg find -x $1 -p hereford.28:25000000-35000000 -c 3 | ${vg}/vg view -gp - > ${bname}.subset.gfa
${vg}/vg find -x $1 -p hereford.28:25000000-35000000 -c 3 | ${vg}/vg view -dp - | ${gviz}/dot -Tsvg -o subgraph1.svg -



