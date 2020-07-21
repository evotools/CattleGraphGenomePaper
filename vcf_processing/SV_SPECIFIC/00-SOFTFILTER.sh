#!/bin/bash

##Selection of large Variants

#VCFs for samples of N'Dama breed were softfiltered then combined using bcftools merge.

for i in Angus01 Angus32065 Angus34122 NDama_ND21 NDama_ND23 NDama_ND39 Sahiwal_3 Sahiwal_6 Sahiwal_7; do
	if [ ! -e RAW/${i}/${i}.vcf.gz.tbi ]; then tabix -p vcf RAW/${i}/${i}.vcf.gz; fi
	bcftools view -O z -f '.,PASS' RAW/${i}/${i}.vcf.gz > RAW/${i}/${i}.filt.vcf.gz && tabix -p vcf RAW/${i}/${i}.filt.vcf.gz && rm RAW/${i}/${i}.vcf.gz*
	echo `pwd`/RAW/${i}/${i}.filt.vcf.gz 
done > samples.txt

