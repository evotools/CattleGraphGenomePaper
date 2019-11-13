#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N graphPrune
#$ -cwd
#$ -l h_rt=23:00:00
#$ -pe sharedmem 4
#$ -R y
#$ -l h_vmem=32.0G
#$ -hold_jid graphAdd,graphIds,graphMap
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
while getopts ":g:n" opt; do
  case $opt in
    g) gv=${OPTARG};;
    n) nversion=${OPTARG};;
  esac
done




export TMPDIR=`pwd`

for i in {1..29}; do
	graph=`head -n $i $gv | tail -1`
	gName=`basename -s ".vg" $graph`
	${vg}/vg prune -t $NSLOTS -u -a -m GRAPH/mapping GRAPH/${gName}.final.vg > GRAPH/chr${i}.final.pruned.vg
	echo "Pruned chr $i"
done
