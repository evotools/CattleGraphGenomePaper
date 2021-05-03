#!/bin/bash

flank=1000
export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/UCSC/
for i in angus ndama ankole brahman; do 
    if [ ! -e GENOMES/${i}.masked.2bit ]; then faToTwoBit GENOMES/${i}.masked.fa GENOMES/${i}.masked.2bit; fi 
    twoBitInfo -nBed GENOMES/${i}.masked.2bit stdout; 
done | \
    awk -v var=$flank 'BEGIN{OFS="\t"}; $2-var < 0{print $1,"0",$3+var}; $2-var >= 0{print $1,$2-var,$3+var}' | \
    bedSort stdin ./GAPS/gaps.bed

