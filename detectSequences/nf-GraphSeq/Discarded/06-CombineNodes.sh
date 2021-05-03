#!/bin/bash

intervals=$1
oname=`basename -s '.bed' ${intervals}`

# Example of header
# Add region lengths
awk 'BEGIN{OFS="\t"};{print $0, $3-$2}' $intervals > node_analysis/${oname}.lengths.bed
# Create output header
echo | awk 'BEGIN{OFS="\t"};{print "SEQID","BPI","BPE","NODES","N_NODES","STRANDS","SEQS","N_CLOSE_TO_GAPS","NODES_LENGTH"}' > node_analysis/${oname}.lengths.merged.bed
# Perform region merging
bedtools merge -c 4,4,5,6,8,9 -o collapse,count,collapse,collapse,sum,sum -d 5 -i node_analysis/${oname}.lengths.bed >> node_analysis/${oname}.lengths.merged.bed
