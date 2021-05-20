#!/bin/bash

# Extract annotate VCF variant type and extract big one only
vcfbreakmulti JOINT/JOINED.SV500.vcf.gz | \
    vcf-annotate --fill-type | \
    python SVLEN.py - 500 > JOINT/JOINED.SV500.annotated.vcf

# Get variant types
if [ ! -e JOINT/CLASSES ] ; then mkdir JOINT/CLASSES; fi
for svtype in `java -jar ~/Documents/Software/snpEff/snpEff/SnpSift.jar extractFields JOINT/JOINED.SV500.annotated.vcf CHROM POS ID TYPE | awk 'NR>1 {print $NF}' | sort | uniq`; do
    fldr=`echo $svtype | tr a-z A-Z`
    if [ ! -e JOINT/CLASSES/$fldr ]; then mkdir JOINT/CLASSES/$fldr; fi
    java -jar ~/Documents/Software/snpEff/snpEff/SnpSift.jar filter "TYPE == '$svtype'" JOINT/JOINED.SV500.annotated.vcf | \
        vcfcreatemulti | \
        java -jar ~/Documents/Software/snpEff/snpEff/SnpSift.jar extractFields - CHROM POS REF ID | \
        python SCRIPTS/ChooseSize.py > JOINT/CLASSES/$fldr/${fldr}.bed 
done