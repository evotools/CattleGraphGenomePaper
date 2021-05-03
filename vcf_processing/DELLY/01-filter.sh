#!/bin/bash

for i in Angus01 Angus32065 Angus34122 ND21 ND23 ND39 Sahiwal_3 Sahiwal_6 Sahiwal_7 ; do 
    bcftools view -f ".,PASS" ${i}/${i}.bqsr.sv.bcf | \
        bcftools filter -e INFO/IMPRECISE=1 -O v > ${i}/${i}.bqsr.sv.filter.precise.vcf
    echo ${i}/${i}.bqsr.sv.filter.precise.vcf
done > sample.filt.txt