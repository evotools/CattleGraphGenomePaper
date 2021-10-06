#!/bin/bash
#$ -cwd 
#$ -N rank
#$ -l h_vmem=48G
#$ -R y
#$ -r y
#$ -l h_rt=23:00:00
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
. /etc/profile.d/modules.sh
#module load roslin/singularity
module load roslin/gcc/7.3.0
module load roslin/python/3.6.8
source ./FORGe-1.1.1/py3/bin/activate

forge_path=${PWD}/FORGe-1.1.1/

method=$1
ref=$2

if [ ! -e RANKED ]; then mkdir RANKED ; fi
if [ ! -e RANKED/$method ]; then mkdir RANKED/$method; fi
${forge_path}/src/rank.py \
	 --method $method --reference $ref.${SGE_TASK_ID}.fa --vars 1KSNP/variants.${SGE_TASK_ID}.1ksnp --window-size 100 --prune 15 --output RANKED/$method/ordered_${method}.${SGE_TASK_ID}.txt
