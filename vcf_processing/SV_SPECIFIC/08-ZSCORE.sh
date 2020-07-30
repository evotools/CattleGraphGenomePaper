#!/bin/bash

for i in Angus NDama Sahiwal; do
    gunzip -c ${i}_final/${i}.consensus.VEP.QUAL30.20DP90.AC5.SV500.vcf.gz | \
        vcfcreatemulti |\
        java -jar ~/Documents/Software/snpEff/snpEff/SnpSift.jar extractFields - CHROM POS REF ID | \
        python SCRIPTS/ChooseSize.py > INTERSECT_BN_ASM/${i}/${i}.consensus.VEP.QUAL30.20DP90.AC5.SV500.bed
    awk '{print $3-$2}' INTERSECT_BN_ASM/${i}/${i}.consensus.VEP.QUAL30.20DP90.AC5.SV500.bed > INTERSECT_BN_ASM/${i}/lengths.txt

    bedtools intersect -wa -a INTERSECT_BN_ASM/${i}/${i}.consensus.VEP.QUAL30.20DP90.AC5.SV500.bed -b INTERSECT_BN_ASM/OMregions.bed | \
        bedtools sort -i - | uniq | bgzip -c > INTERSECT_BN_ASM/${i}/${i}.consensus.VEP.QUAL30.20DP90.AC5.SV500.intersectOM.bed.gz && \
        tabix -f -p bed INTERSECT_BN_ASM/${i}/${i}.consensus.VEP.QUAL30.20DP90.AC5.SV500.intersectOM.vcf.gz
    bedtools subtract -A -a INTERSECT_BN_ASM/${i}/${i}.consensus.VEP.QUAL30.20DP90.AC5.SV500.bed -b INTERSECT_BN_ASM/OMregions.bed | \
        bedtools sort -i - | uniq | bgzip -c > INTERSECT_BN_ASM/${i}/${i}.consensus.VEP.QUAL30.20DP90.AC5.SV500.notInOM.bed.gz && \
        tabix -f -p bed INTERSECT_BN_ASM/${i}/${i}.consensus.VEP.QUAL30.20DP90.AC5.SV500.notInOM.vcf.gz

    CWD=$PWD
    cd INTERSECT_BN_ASM/${i}
    cp $CWD/INTERSECT_BN_ASM_FULLSVset/ARS_UCD1.2.fai ./
    nconfirmed=`gunzip -c ${i}.consensus.VEP.QUAL30.20DP90.AC5.SV500.intersectOM.bed.gz | awk 'BEGIN{n=0}; {n+=1}; END{print n}'`
    cp $CWD/INTERSECT_BN_ASM/OMregions.bed ./
    Rscript ${CWD}/INTERSECT_BN_ASM_FULLSVset/TestRandomRegions.R ${nconfirmed} 10000
    cd ../../
done