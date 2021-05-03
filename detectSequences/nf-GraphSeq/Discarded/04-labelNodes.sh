#!/bin/bash
intervals=$1
nmers=$2
bname=`basename -s '.bed' $intervals`

#bedtools intersect -a $intervals -b $nmers -wa | awk '{print $4}' > intersected.txt
bedtools intersect -a $intervals -b $nmers -v | awk 'BEGIN{OFS="\t"};{print $0, "0"}' > ${bname}.noOverlaps.bed
bedtools intersect -a $intervals -b $nmers -u | awk 'BEGIN{OFS="\t"};{print $0, "1"}' > ${bname}.Overlaps.bed
cat ${bname}.noOverlaps.bed ${bname}.Overlaps.bed | bedSort stdin ${bname}.labeled.bed
