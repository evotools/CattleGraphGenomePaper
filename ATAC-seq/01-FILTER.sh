#!/bin/bash
# Script to process peaks in each file

for w in HEREFORD ANGUS NDAMA ANKOLE BRAHMAN PANGENOME; do
	genome_len=`grep -w $w GENOME_SIZES.txt | awk '{print $2}'`;
	genome_name=`python -c "import sys; print(sys.argv[1].lower())" $w `

	if [ ! -e ${w}/BLANKPEAKS/ ]; then mkdir ${w}/BLANKPEAKS/; fi

	for blank in HF3457_nucfree_ATAC; do
		if [ ! -e ${w}/BLANKPEAKS/${blank} ]; then
	                mkdir ${w}/BLANKPEAKS/${blank};
	        fi;
	        awk 'BEGIN{OFS="\t"};{print $1,$2,$3,$4"#"$9}' ${w}/GENRICH/${blank}/*.narrowPeak > ${w}/BLANKPEAKS/${blank}/BLANKPEAKS.bed
		echo "Blank ready."
	done

	for i in HF3457_Bcell_ATAC ND230_Bcell_ATAC Nelore1_Bcell_ATAC; do 
		bedtools subtract -A -f 0.50 -a ${w}/GENRICH/${i}/${i}.narrowPeak -b ${w}/BLANKPEAKS/${blank}/BLANKPEAKS.bed | sort | uniq > ${w}/GENRICH/${i}/${i}.noblank.narrowPeak
		bedtools subtract -A -f 0.50 -a ${w}/GENRICH/${i}/${i}.noblank.narrowPeak -b ./SOFTMASKED/${genome_name}.softmasked.bed | sort | uniq > ${w}/GENRICH/${i}/${i}.noblank.lowrep.narrowPeak
		Rscript SCRIPTS/QSCORES.R ${w}/GENRICH/${i}/${i}.noblank.lowrep.narrowPeak $genome_len
		Rscript SCRIPTS/QSCORES.R ${w}/GENRICH/${i}/${i}.noblank.narrowPeak $genome_len
		Rscript SCRIPTS/QSCORES.R ${w}/GENRICH/${i}/${i}.narrowPeak $genome_len
		echo "Initial peaks                 : "`wc -l ${w}/GENRICH/${i}/${i}.narrowPeak`
		echo "No-blank peaks                : "`wc -l ${w}/GENRICH/${i}/${i}.noblank.narrowPeak`
		echo "Low repetitive elements peaks : "`wc -l ${w}/GENRICH/${i}/${i}.noblank.lowrep.narrowPeak`
	done
	echo "Done ${w}"
done
