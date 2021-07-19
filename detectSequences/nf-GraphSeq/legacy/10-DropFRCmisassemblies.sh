#!/bin/bash
ctgs=$1
frc_feat=$2
outname=$3
# Example input node_analysis/candidate_75/nodes_candidate_75.bed
# Example of outname node_analysis/candidate_50/nodes_candidate_75
bedtools intersect -a ${frc_feat} -b ${ctgs} -wa -wb | \
    awk 'BEGIN{OFS="\t"}; $3=="HIGH_OUTIE_PE" || $3=="HIGH_SINGLE_PE" || $3=="HIGH_SPAN_PE" || $3~"LOW_" {print $10, $11,$12, $12-$11}' | \
    sort | \
    uniq > ${outname}.remove_by_FRC.bed
bedtools intersect -a ${ctgs} -b ${outname}.remove_by_FRC.bed -v > ${outname}.post_FRC.bed