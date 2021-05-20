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

if [ ! -e BNOM/BED_V2 ]; then mkdir BNOM/BED_V2; else rm BNOM/BED_V2/*/*.bed; fi

for sample in NDama_1_1628 NDama_NN031; do
    # Create folder structure if it doesn't exists
    for svtype in `awk '$1!~"#"{print $5}' BNOM/VCFs/${sample}_fix.vcf | sort | uniq | sed 's/<//g' | sed 's/>//g'`; do 
        if [ ! -e BNOM/BED_V2/${svtype} ]; then mkdir BNOM/BED_V2/${svtype}; fi;  

        # Extract SVs for class
        bcftools query -f'%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\t%CIPOS\t%CIEND\n' BNOM/VCFs/${sample}_fix.vcf.gz | \
            awk -v val=${svtype} '$5==val || $1~"#" {print}' | \
            python SCRIPTS/ExpandRegion.py 1000 | \
            sort -k1,1n -k2,2n -k3,3n >> BNOM/BED_V2/${svtype}/${svtype}.bed    
    done
done

for foldername in `ls BNOM/BED_V2/`; do 
    bname=`basename $foldername`
    if [ $bname != "TRA" ]; then
        bedtools sort -i BNOM/BED_V2/${bname}/${bname}.bed | bedtools merge -i - > BNOM/BED_V2/${bname}/${bname}.sort.bed
        echo $bname
    fi
done