#!/bin/bash

gunzip -c JOINT/JOINED.SV500.vcf.gz | \
    vcfcreatemulti | \
    java -jar ~/Documents/Software/snpEff/snpEff/SnpSift.jar extractFields - CHROM POS REF ID | \
    python SCRIPTS/ChooseSize.py > INTERSECT_BN_ASM_FULLSVset/JOINED.SV500.bed
awk '{print $3-$2}' INTERSECT_BN_ASM_FULLSVset/JOINED.SV500.bed > INTERSECT_BN_ASM_FULLSVset/lengths.txt

bedtools intersect -wa -a INTERSECT_BN_ASM_FULLSVset/JOINED.SV500.bed -b INTERSECT_BN_ASM/OMregions.bed | \
    bedtools sort -i - | uniq | bgzip -c > INTERSECT_BN_ASM_FULLSVset/JOINED.SV500.intersectOM.bed.gz && \
    tabix -f -p bed INTERSECT_BN_ASM_FULLSVset/JOINED.SV500.intersectOM.bed.gz
bedtools subtract -A -a INTERSECT_BN_ASM_FULLSVset/JOINED.SV500.bed -b INTERSECT_BN_ASM/OMregions.bed | \
    bedtools sort -i - | uniq | bgzip -c > INTERSECT_BN_ASM_FULLSVset/JOINED.SV500.notInOM.bed.gz && \
    tabix -f -p bed INTERSECT_BN_ASM_FULLSVset/JOINED.SV500.notInOM.bed.gz

CWD=$PWD
cd INTERSECT_BN_ASM_FULLSVset
nconfirmed=`gunzip -c ${CWD}/INTERSECT_BN_ASM_FULLSVset/JOINED.SV500.intersectOM.bed.gz | awk 'BEGIN{n=0}; {n+=1}; END{print n}'`
cp $CWD/INTERSECT_BN_ASM/OMregions.bed ./

# Run with nextflow to keep a high throughput on eddie
nextflow run main.nf --faifile $PWD/ARS_UCD1.2.fai --ntest 10000 --omregions $PWD/OMregions.bed --lengths $PWD/lengths.txt --positives $nconfirmed
cd ..

