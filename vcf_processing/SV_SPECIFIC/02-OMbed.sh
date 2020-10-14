#!/bin/bash

if [ ! -e BNOM ]; then mkdir BNOM; fi
if [ ! -e BNOM/BED ]; then mkdir BNOM/BED; fi


for sample in NDama_1_1628 NDama_NN031; do
    bcftools query -f'%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\t%CIPOS\t%CIEND\n' BNOM/VCFs/${sample}_fix.vcf.gz | \
        grep -v "TRA" | 
        python SCRIPTS/ExpandRegion.py 1000 | \
        sort -k1,1n -k2,2n -k3,3n > BNOM/BED/${sample}.SV.noTranslocations.bed
done

# Extract the ones into a region from OM
cat $( for sample in NDama_1_1628 NDama_NN031; do echo BNOM/BED/${sample}.SV.noTranslocations.bed; done ) | \
    bedtools sort -i - | \
    bedtools merge -i - | \
    awk '$1>=1 && $1<=29 {print}' > INTERSECT_BN_ASM/OMregions.bed
