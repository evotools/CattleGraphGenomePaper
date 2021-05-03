#!/bin/bash

# Add vcflib
export PATH=$PATH:~/Documents/Software/vcflib/bin

## Repeat for Angus samples
# Get Angus specific variants removing the ones found in other samples 
for s in Angus01 Angus32065 Angus34122; do
	input=RAW/${s}/${s}.filt.vcf.gz
	bname=${s}
	for i in NDama_ND21 NDama_ND23 NDama_ND39 Sahiwal_3 Sahiwal_6 Sahiwal_7; do
		bcftools isec -p Angus_isec -O z $input RAW/${i}/${i}.filt.vcf.gz && 	rm  ./${s}_*.vcf.gz*
		bname=${bname}_${i}
		mv Angus_isec/0000.vcf.gz ./${bname}.vcf.gz && tabix -p vcf ${bname}.vcf.gz && rm -rf Angus_isec/
		input=${bname}.vcf.gz
		echo "Filtered $s with $i"
	done
	if [ ! -e Angus_final ]; then mkdir Angus_final/ ; fi
	if [ ! -e Angus_final/$s ]; then mkdir Angus_final/$s ; fi
	mv ${input}  Angus_final/${s}/${s}.unique.vcf.gz && tabix -p vcf Angus_final/${s}/${s}.unique.vcf.gz
	echo `pwd`/Angus_final/${s}/${s}.unique.vcf.gz >> Angus.final.txt
done
# Combine filtered N'Dama
bcftools merge -m all -l Angus.final.txt -O v | sed 's/hereford.//g' | python SVLEN.py - | bgzip -c > Angus_final/Angus.vcf.gz && tabix -p vcf Angus_final/Angus.vcf.gz

## Annotate variants
mkdir VEP/output 
cp Angus_final/Angus.vcf.gz ${HOME}/vep_data/input
docker run -t -i -v $HOME/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep vep -i /opt/vep/.vep/input/Angus.vcf.gz -o /opt/vep/.vep/output/Angus.consensus.VEP.vcf.gz \
	--vcf --fork 4 --species bos_taurus --variant_class --sift b --nearest symbol --distance 200 --cache --offline --compress_output bgzip \
	--dir_cache /opt/vep/.vep/ \
	--dir_plugins /opt/vep/.vep/Plugins/ > VEP/angus.ensemblv100_ARS-UCD1.2.log
mv ${HOME}/vep_data/output/Angus.vcf_*.html ./VEP/output
mv ${HOME}/vep_data/output/Angus.consensus.VEP.vcf.gz ./VEP/output/ && tabix -p vcf ./VEP/output/Angus.consensus.VEP.vcf.gz 

# Filter N'Dama-specific annotated variants 
bcftools view -f 'QUAL>30' -O v ./VEP/output/Angus.consensus.VEP.vcf.gz | \
	vcftools --vcf - --min-meanDP 20 --max-meanDP 90 --non-ref-ac 5 --recode --recode-INFO-all --stdout |
	python SVLEN.py - 500 |
	bgzip -c > Angus_final/Angus.consensus.VEP.QUAL30.20DP90.AC5.SV500.vcf.gz && \
	tabix -p vcf Angus_final/Angus.consensus.VEP.QUAL30.20DP90.AC5.SV500.vcf.gz
	
# Intersect with OMdata
if [ ! -e INTERSECT_BN_ASM/Angus ]; then mkdir INTERSECT_BN_ASM/Angus ; fi
bedtools intersect -header -a Angus_final/Angus.consensus.VEP.QUAL30.20DP90.AC5.SV500.vcf.gz -b INTERSECT_BN_ASM/OMregions.bed | vcf-sort | vcfuniq | bgzip -c > INTERSECT_BN_ASM/Angus/Angus.consensus.VEP.QUAL30.20DP90.AC5.SV500.intersectOM.vcf.gz
bedtools subtract -header -a Angus_final/Angus.consensus.VEP.QUAL30.20DP90.AC5.SV500.vcf.gz -b INTERSECT_BN_ASM/OMregions.bed | vcf-sort | vcfuniq |bgzip -c > INTERSECT_BN_ASM/Angus/Angus.consensus.VEP.QUAL30.20DP90.AC5.SV500.notInOM.vcf.gz

echo "Done Angus"