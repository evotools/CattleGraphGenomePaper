#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N graphIndex
#$ -cwd
#$ -l h_rt=8:00:00
#$ -pe sharedmem 4
#$ -R y
#$ -l h_vmem=8.0G
#$ -hold_jid gc_60419
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
graph=`head -n $SGE_TASK_ID $1 | tail -1`
vcf=`head -n $SGE_TASK_ID $2 | tail -1`
spp=$3
gName=`basename -s ".vg" $1`
export TMPDIR=`pwd`

cd GRAPH
${vg}/vg ids -j $(for i in $(seq 1 29); do echo ${gName}.final.vg; done)
${vg}/vg index -p -x GRAPH/CHR${SGE_TASK_ID}/${gName}.final.xg GRAPH/CHR${SGE_TASK_ID}/${gName}.final.vg
${vg}/vg prune -p -r GRAPH/CHR${SGE_TASK_ID}/${gName}.final.vg > GRAPH/CHR${SGE_TASK_ID}/${gName}.final.prune.vg
${vg}/vg index -p -g GRAPH/CHR${SGE_TASK_ID}/${gName}.final.gcsa GRAPH/CHR${SGE_TASK_ID}/${gName}.final.prune.vg
echo "Done chr ${SGE_TASK_ID}"

