#!/bin/bash
if [ ! -e VENN ]; then mkdir VENN; fi
awk '{print "hereford."$0}' INTERSECT_BN_ASM_FULLSVset/JOINED.SV500.bed > VENN/VG5P.bed
awk '{print "hereford."$0}' INTERSECT_BN_ASM_FULLSVset/OMregions.bed > VENN/BIONANO.bed
cp SV_CALLERS/MERGED/ALL/delly.merged.min500bp.1suppcall.100bp_breakpoint.noTRA.isecBNOM.vcf > VENN/DELLY.bed

cd VENN
if [ ! -e RES ]; then mkdir RES; fi
# Specific
bedtools intersect -header -a BIONANO.bed -b DELLY.vcf -v -header | bedtools intersect -a - -b VG5P.bed -v > RES/uniqueBN.bed
bedtools intersect -header -a DELLY.vcf -b BIONANO.bed -v -header | bedtools intersect -a - -b VG5P.bed -v > RES/uniqueDE.bed
bedtools intersect -header -a VG5P.bed -b BIONANO.bed -v -header | bedtools intersect -a - -b DELLY.vcf -v > RES/uniqueVG.bed
echo BIONANO `wc -l RES/uniqueBN.bed`
echo DELLY `wc -l RES/uniqueDE.bed`
echo VG5P `wc -l RES/uniqueVG.bed`

# shared in pairs
bedtools intersect -header -a BIONANO.bed -b DELLY.vcf -header -u | bedtools intersect -a - -b VG5P.bed -v > RES/BNxDE.bed
echo BIONANO DELLY `wc -l RES/BNxDE.bed`
bedtools intersect -header -a BIONANO.bed -b VG5P.bed -header -u | bedtools intersect -a - -b DELLY.vcf -v > RES/BNxVG.bed
echo BIONANO VG5P `wc -l RES/BNxVG.bed`
bedtools intersect -header -a DELLY.vcf -b VG5P.bed -header -u | bedtools intersect -a - -b BIONANO.bed -v > RES/DExVG.bed
echo DELLY VG5P `wc -l RES/DExVG.bed`
bedtools intersect -header -a DELLY.vcf -b BIONANO.bed -header -u | bedtools intersect -a - -b VG5P.bed -v > RES/DExBN.bed
echo DELLY BIONANO `wc -l RES/DExBN.bed`
bedtools intersect -header -a VG5P.bed -b BIONANO.bed -header -u | bedtools intersect -a - -b DELLY.vcf -v > RES/VGxBN.bed
echo VG5P BIONANO `wc -l RES/VGxBN.bed`
bedtools intersect -header -a VG5P.bed -b DELLY.vcf -header -u | bedtools intersect -a - -b BIONANO.bed -v > RES/VGxDE.bed
echo VG5P DELLY `wc -l RES/VGxDE.bed`

# All three
bedtools intersect -a BIONANO.bed -b DELLY.vcf -header -u | bedtools intersect -a - -b VG5P.bed -u -header > RES/BNall.bed
echo BIONANO DELLY VG5P `wc -l RES/BNall.bed`
bedtools intersect -a DELLY.vcf -b BIONANO.bed -header -u | bedtools intersect -a - -b VG5P.bed -u -header > RES/DEall.bed
echo DELLY BIONANO VG5P `wc -l RES/DEall.bed`
bedtools intersect -a VG5P.bed -b BIONANO.bed -header -u | bedtools intersect -a - -b DELLY.vcf -u -header > RES/VGall.bed
echo VG5P DELLY BIONANO `wc -l RES/VGall.bed`

