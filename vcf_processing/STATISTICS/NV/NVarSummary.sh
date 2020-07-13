#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N ABsum
#$ -cwd
#$ -l h_rt=00:20:00
#$ -R y
#$ -l h_vmem=8.0G
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh
. /etc/profile.d/modules.sh
# module load roslin/bcftools
# module load igmm/apps/tabix
module load anaconda
module load roslin/gcc
source activate DataPy
module load igmm/apps/tabix

sample=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $1}'`
vcf=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $2}'`
isectgt=`realpath $2`
tgtfld=$3
outfld=$4
currpath=`pwd`


python ${currpath}/NormaliseQUAL.py -i ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.tmp.gz | bgzip -c > ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.txt.gz && rm ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.tmp.gz
python ${currpath}/SummariseNV.py ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.txt.gz > ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.summary.csv

python ${currpath}/NormaliseQUAL.py -i ${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.tmp.gz | bgzip -c > ${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.txt.gz && rm ${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.tmp.gz
python ${currpath}/SummariseNV.py ${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.txt.gz > ${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.summary.csv

python ${currpath}/NormaliseQUAL.py -i ${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.tmp.gz | bgzip -c > ${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.txt.gz && rm ${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.tmp.gz
python ${currpath}/SummariseNV.py ${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.txt.gz > ${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.summary.csv
