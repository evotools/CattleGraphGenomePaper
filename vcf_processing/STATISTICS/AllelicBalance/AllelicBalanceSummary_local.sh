#!/bin/bash
module load roslin/bcftools
module load igmm/apps/tabix
module load anaconda
source activate DataPy

while read p; do 
sample=`echo ${p} | awk '{print $1}'`
vcf=`echo ${p} | awk '{print $2}'`
isectgt=`realpath $2`
tgtfld=$3
outfld=$4
currpath=`pwd`


python ${currpath}/SummariseAB.py ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.txt.gz > ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.summary.csv
python ${currpath}/SummariseAB.py ${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.txt.gz > ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.summary.csv
python ${currpath}/SummariseAB.py ${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.txt.gz > ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.summary.csv

done < $1
