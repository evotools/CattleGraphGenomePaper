#!/bin/bash

for i in ANGUS NDAMA SAHIWAL; do 
    awk '$1~"#"{print};$1!~"#" && $0!~"SVTYPE=TRA"{print}' ../MERGED/$i/$i.SV500.vcf | \
        bedtools intersect -a - -b ../regions.bed -wa -u -header  > ../MERGED/$i/$i.SV500.noTRA.isecBNOM.vcf 
done
awk '$1~"#"{print};$1!~"#" && $0!~"SVTYPE=TRA"{print}' ../MERGED/ALL/delly.merged.min500bp.1suppcall.100bp_breakpoint.vcf |
    bedtools intersect -a - -b ../regions.bed -wa -u -header > ../MERGED/ALL/delly.merged.min500bp.1suppcall.100bp_breakpoint.noTRA.isecBNOM.vcf