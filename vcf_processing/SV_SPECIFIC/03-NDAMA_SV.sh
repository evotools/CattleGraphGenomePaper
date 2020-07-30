#!/bin/bash

# Add vcflib
export PATH=$PATH:~/Documents/Software/vcflib/bin

## Combine N'Dama samples
# Get N'Dama specific variants removing the ones found in other samples 
for n in NDama_ND21 NDama_ND23 NDama_ND39; do
	input=RAW/${n}/${n}.filt.vcf.gz
	# # Convert to multi-line per allele
	# vcfmultibreak RAW/${n}/${n}.filt.vcf.gz | bgzip -c > RAW/${n}/${n}.filt.multib.vcf.gz && tabix -p vcf RAW/${n}/${n}.filt.multib.vcf.gz
	bname=${n}
	for i in Angus01 Angus32065 Angus34122 Sahiwal_3 Sahiwal_6 Sahiwal_7; do
		bcftools isec -p NDama_isec -O z $input RAW/${i}/${i}.filt.vcf.gz && 	rm  ./${n}_*.vcf.gz*
		# bcftools isec -p NDama_isec -O z $input RAW/${i}/${i}.filt.multib.vcf.gz && 	rm  ./${n}_*.vcf.gz*
		bname=${bname}_${i}
		mv NDama_isec/0000.vcf.gz ./${bname}.vcf.gz && tabix -p vcf ${bname}.vcf.gz && rm -rf NDama_isec/
		input=${bname}.vcf.gz
		echo "Filtered $n with $i"
	done && rm RAW/${n}/${n}.filt.multib.vcf.gz RAW/${n}/${n}.filt.multib.vcf.gz.tbi
	if [ ! -e NDama_final ]; then mkdir NDama_final/ ; fi
	if [ ! -e NDama_final/$n ]; then mkdir NDama_final/$n ; fi
	mv ${input}  NDama_final/${n}/${n}.unique.vcf.gz && tabix -p vcf NDama_final/${n}/${n}.unique.vcf.gz
	echo `pwd`/NDama_final/${n}/${n}.unique.vcf.gz > Ndama.final.txt
done
# Combine filtered N'Dama
bcftools merge -l NDama.final.txt -O v | sed 's/hereford.//g' | python SVLEN.py - | bgzip -c > NDama_final/NDama.vcf.gz && tabix -p vcf NDama_final/NDama.vcf.gz
# bcftools merge -l NDama.final.txt -O v | vcfmultibreak | sed 's/hereford.//g' | python SVLEN.py - | bgzip -c > NDama_final/NDama.vcf.gz && tabix -p vcf NDama_final/NDama.vcf.gz

## Annotate variants
mkdir VEP/output 
ln -s NDama_final/NDama.vcf.gz ${HOME}/vep_data/input
docker run -t -i -v $HOME/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep vep -i /opt/vep/.vep/input/NDama.vcf.gz -o /opt/vep/.vep/output/NDama.vcf \
	--vcf --fork 4 --species bos_taurus --variant_class --sift b --nearest symbol --distance 200 --cache --offline \
	--dir_cache /opt/vep/.vep/ \
	--dir_plugins /opt/vep/.vep/Plugins/ > VEP/ensemblv100_ARS-UCD1.2.log
mv ${HOME}/vep_data/output/NDama.vcf_*.html ./VEP/output
cat ${HOME}/vep_data/output/NDama.vcf | bgzip -c > ./VEP/output/NDama.consensus.VEP.vcf.gz && tabix -p vcf ./VEP/output/NDama.consensus.VEP.vcf.gz && rm ${HOME}/vep_data/output/* ${HOME}/vep_data/input/*

# Filter N'Dama-specific annotated variants 
bcftools view -f 'QUAL>30' -O v ./VEP/output/NDama.consensus.VEP.vcf.gz | \
	vcftools --vcf - --min-meanDP 20 --max-meanDP 90 --non-ref-ac 5 --recode --recode-INFO-all --stdout |
	python SVLEN.py - 500 |
	bgzip -c > NDama_final/NDama.consensus.VEP.QUAL30.20DP90.AC5.SV500.vcf.gz && \
	tabix -p vcf NDama_final/NDama.consensus.VEP.QUAL30.20DP90.AC5.SV500.vcf.gz

# Intersect with OMdata
if [ ! -e INTERSECT_BN_ASM/NDama ]; then mkdir INTERSECT_BN_ASM/NDama ; fi
bedtools intersect -header -a NDama_final/NDama.consensus.VEP.QUAL30.20DP90.AC5.SV500.vcf.gz -b INTERSECT_BN_ASM/OMregions.bed | vcf-sort | vcfuniq | bgzip -c > INTERSECT_BN_ASM/NDama.consensus.VEP.QUAL30.20DP90.AC5.SV500.intersectOM.vcf.gz
bedtools subtract -header -a NDama_final/NDama.consensus.VEP.QUAL30.20DP90.AC5.SV500.vcf.gz -b INTERSECT_BN_ASM/OMregions.bed | vcf-sort | vcfuniq | bgzip -c > INTERSECT_BN_ASM/NDama.consensus.VEP.QUAL30.20DP90.AC5.SV500.notInOM.vcf.gz

echo "Done N'Dama"