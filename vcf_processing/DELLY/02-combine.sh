#!/bin/bash

head -3 sample.filt.txt > angus.txt
head -6 sample.filt.txt | tail -3 > ndama.txt
tail -3 sample.filt.txt > sahiwal.txt

for i in angus.txt ndama.txt sahiwal.txt; do 
    bname=`echo $i | sed 's/.txt//g' | awk '{print toupper($1)}'`
    if [ ! -e ../MERGED/${bname} ]; then mkdir ../MERGED/${bname}; fi 
    ~/Documents/Software/SURVIVOR/Debug/SURVIVOR merge ${i} 100 1 1 0 0 500 ../MERGED/${bname}/${bname}.SV500.vcf
    vcftools --vcf ../MERGED/${bname}/${bname}.SV500.vcf --max-missing 1 --recode --recode-INFO-all --stdout | \
        bcftools view -c 5:nref  > ../MERGED/${bname}/${bname}.SV500.minAC5.nomiss.vcf
done

if [ ! -e ../MERGED/ALL ]; then mkdir ../MERGED/ALL; fi 
~/Documents/Software/SURVIVOR/Debug/SURVIVOR merge sample.filt.txt 100 1 1 0 0 500 ../MERGED/ALL/delly.merged.min500bp.1suppcall.100bp_breakpoint.vcf