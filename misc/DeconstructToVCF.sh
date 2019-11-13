#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N graphDec
#$ -cwd
#$ -l h_rt=47:00:00
#$ -pe sharedmem 1
#$ -R y
#$ -l h_vmem=256.0G
#$ -hold_jid gc_60419
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.err
#$ -P roslin_ctlgh
. /etc/profile.d/modules.sh
module load java/jdk/1.8.0
module load roslin/samtools/1.9
module load igmm/apps/tabix

vg=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/vg
gviz=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/GVIZ/bin/

# 
# Start program
# 

if [ ! -e all.snarls ]; then
	$vg/vg snarls GRAPH/all.xg > all.snarls
fi
$vg/vg deconstruct -p hereford.28 -r all.snarls -t $NSLOTS GRAPH/all.xg | bgzip -c > hereford.28.vcf.gz


