#!/bin/bash
intervals=$1
CE=$2
repeatMasker=$3

# Add details about CE
bedtools intersect -a ${intervals} -b ${CE} -wa -wb > withCE.bed

# Add details about repetitive elements
bedtools intersect -a ${intervals} -b ${repeatMasker} -wa -wb > withRM.bed


