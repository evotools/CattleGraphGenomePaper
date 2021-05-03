#!/bin/bash

intervals=$1
outf=$2
name=$3
threshold=$4

if [ ! -e $outf ]; then mkdir $outf; fi
awk 'NR==1 {print}; NR>1 && $11~"LONG"{print}' ${intervals} > ${outf}/${name}.long.bed
awk -v val=$threshold 'NR==1{print};NR>1 && $9/$10>val {print}' ${outf}/${name}.long.bed > ${outf}/${name}.long.novel.bed
awk 'NR==1{print};NR>1 && $11!~"TELOMER"{print}' ${outf}/${name}.long.novel.bed > ${outf}/${name}.long.novel.noTelomere.bed
awk 'NR==1{print};NR>1 && $11!~"FLANK"{print}' ${outf}/${name}.long.novel.noTelomere.bed > ${outf}/${name}.long.novel.noTelomere.noFlankGaps.bed