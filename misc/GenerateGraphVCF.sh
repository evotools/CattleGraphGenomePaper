#!/bin/bash
#$ -cwd
#$ -N graphVCF
#$ -pe sharedmem 1
#$ -R y
#$ -l h_rt=47:59:59
#$ -l h_vmem=16G
#$ -P roslin_ctlgh
#$ -e LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -o LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out

## Script to generate a vg compliant graph from a vcf file


. /etc/profile.d/modules.sh
module load anaconda
module load roslin/bedtools
module load igmm/apps/vcftools
module load igmm/apps/tabix
module load R
source activate DataPy3


vcftools --gzvcf $1 --chr $SGE_TASK_ID --recode --recode-INFO-all --stdout | python ./GraphVCF.py - stdout | awk '$1~"#"{print}; $1!~"#" {print "hereford."$0}'| bgzip -c > GRAPH_${SGE_TASK_ID}.vcf.gz
tabix -p vcf GRAPH_${SGE_TASK_ID}.vcf.gz

