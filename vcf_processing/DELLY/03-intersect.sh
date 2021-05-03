#!/bin/bash

for i in ANGUS NDAMA SAHIWAL; do 
    bedtools intersect -a ../MERGED/$i/$i.SV500.vcf -b ../regions.bed -wa -u -header > ../MERGED/$i/$i.SV500.isecBNOM.vcf 
done

bedtools intersect -a ../MERGED/ALL/delly.merged.min500bp.1suppcall.100bp_breakpoint.vcf -b ../regions.bed -wa -u -header > ../MERGED/ALL/delly.merged.min500bp.1suppcall.100bp_breakpoint.isecBNOM.vcf