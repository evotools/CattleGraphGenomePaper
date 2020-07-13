#!/bin/bash
# Script to process peaks in each file

for w in HEREFORD ANGUS NDAMA ANKOLE BRAHMAN PANGENOME; do
	genome_len=`grep -w $w GENOME_SIZES.txt | awk '{print $2}'`;

	if [ ! -e ${w}/NONHEREFORD/ ]; then mkdir ${w}/NONHEREFORD/; fi

	for i in HF3457_nucfree_ATAC; do
		if [ ! -e ${w}/NONHEREFORD/${i} ]; then
	                mkdir ${w}/NONHEREFORD/${i};
	        fi;
	        awk 'BEGIN{OFS="\t"};{print $1,$2,$3,$4"#"$9}' ${w}/GENRICH/${i}/*.narrowPeak > ${w}/NONHEREFORD/${i}/NONHEREFORD.bed
		echo "Blank ready."
	done

	for i in HF3457_Bcell_ATAC ND230_Bcell_ATAC Nelore1_Bcell_ATAC; do 
		if [ ! -e ${w}/NONHEREFORD/${i} ]; then 
			mkdir ${w}/NONHEREFORD/${i}; 
		fi; 
		bedtools subtract -A -f 0.50 -a ${w}/GENRICH/${i}/*.narrowPeak -b ${w}/NONHEREFORD/HF3457_nucfree_ATAC/NONHEREFORD.bed | sort | uniq > ${w}/GENRICH/${i}/${i}.noblank.narrowPeak
		Rscript SCRIPTS/QSCORES.R ${w}/GENRICH/${i}/${i}.noblank.narrowPeak $genome_len
		Rscript SCRIPTS/QSCORES.R ${w}/GENRICH/${i}/${i}.narrowPeak $genome_len
	done
	echo "Done ${w}"
done
