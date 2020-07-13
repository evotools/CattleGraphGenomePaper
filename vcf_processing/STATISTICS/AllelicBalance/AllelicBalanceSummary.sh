#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N ABsum
#$ -cwd
#$ -l h_rt=00:20:00
#$ -R y
#$ -l h_vmem=1.0G
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh
. /etc/profile.d/modules.sh
module load anaconda
source activate DataPy

sample=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $1}'`
vcf=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $2}'`
isectgt=`realpath $2`
tgtfld=$3
outfld=$4
currpath=`pwd`


python ${currpath}/SummariseAB.py ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.txt.gz > ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.summary.csv
python ${currpath}/SummariseAB.py ${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.txt.gz > ${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.summary.csv
python ${currpath}/SummariseAB.py ${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.txt.gz > ${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.summary.csv

