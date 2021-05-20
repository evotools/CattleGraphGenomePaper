#!/bin/bash
#$ -N bamcov
#$ -R y
#$ -pe sharedmem 4
#$ -l h_vmem=4G
#$ -l h_rt=23:59:59
#$ -P roslin_ctlgh
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -cwd

. /etc/profile.d/modules.sh
module load anaconda
source activate deeptools

sample=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $1}'`
bamfile=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $2}'`
fbasedir=`dirname $bamfile`
maptype=$2 # UNIQUE/MULTIMAP

samtools view -hb -L ./LISTS/ctg2keep.bed ${bamfile} | \
		samtools sort -T DEEPTOOLS/${maptype}/${sample}/ -@4 -m 10G > DEEPTOOLS/${maptype}/${sample}/${sample}.bam.sort
samtools index -@ 4 DEEPTOOLS/${maptype}/${sample}/${sample}.bam.sort
bamCoverage -b DEEPTOOLS/${maptype}/${sample}/${sample}.bam.sort -o DEEPTOOLS/${maptype}/${sample}/${sample}.bw --minFragmentLength 35 --maxFragmentLength 150 --normalizeUsing RPGC -p $NSLOTS -bs 10 -e --effectiveGenomeSize 2779691414
