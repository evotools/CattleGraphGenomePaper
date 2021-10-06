#!/bin/bash
#$ -cwd 
#$ -N build
#$ -l h_vmem=64G
#$ -R y
#$ -r y
#$ -l h_rt=23:00:00
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
. /etc/profile.d/modules.sh
module load roslin/gcc/7.3.0
module load roslin/python/3.6.8
source ./FORGe/py3/bin/activate

forge_path=${PWD}/FORGe-1.1.1/

method=$1
ref=$2
pct=$3
blowup=$4

# ${forge_path}/src/rank.py --method $method --reference $ref --vars variants.1ksnp --window-size 100 --prune 15 --output RANKED/ordered_${method}.txt
suffix=""
if [ $blowup -eq 1 ]; then suffix=".blowup"; fi

if [ ! -e RES ]; then mkdir RES ; fi
if [ ! -e RES/${method}${suffix} ]; then mkdir RES/${method}${suffix} ; fi
if [ ! -e RES/${method}${suffix}/$pct ]; then mkdir RES/${method}${suffix}/$pct ; fi

${forge_path}/src/build.py \
	--reference ${ref}.${SGE_TASK_ID}.fa \
	--vars ./1KSNP/variants.${SGE_TASK_ID}.1ksnp \
	--window-size 100 \
	--hisat RES/${method}${suffix}/${pct}/variants.${method}${suffix}.${pct}.${SGE_TASK_ID}.snp \
	--sorted RANKED/${method}/ordered_${method}.${SGE_TASK_ID}.txt${suffix} \
	--pct ${pct} 
