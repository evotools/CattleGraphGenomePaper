#!/bin/bash

for i in ANGUS ANKOLE NDAMA BRAHMAN HEREFORD; do 
	for w in `ls ${i}/GENRICH/`; do
		fname=`echo $w | sed 's+/++g'`
		cat ${i}/GENRICH/${w}/*.narrowPeak | awk 'BEGIN{OFS="\t"};{print $1,$2,$3,$4,$5,$6}' > ${i}/GENRICH/${w}/${fname}.bed
		echo "Done $i, $w"
	done
done
