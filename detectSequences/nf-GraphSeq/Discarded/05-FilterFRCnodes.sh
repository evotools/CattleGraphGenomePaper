#!/bin/bash
intervals=$1
frc=$2
bname=`basename -s '.bed' $intervals`

bedtools intersect -a $intervals -b $frc -v | awk 'BEGIN{OFS="\t"};{print $0}' > node_analysis/${bname}.no_frc.bed

# bedtools intersect -a $intervals -b $frc -u | awk '{print $4}' > node_analysis/nodes_removed.txt

