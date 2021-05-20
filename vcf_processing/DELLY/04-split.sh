#!/bin/bash

if [ ! -e SV_CALLERS/MERGED/ALL ]; then mkdir SV_CALLERS/MERGED/ALL; fi
for svtype in `awk '$1!~"#"{print $5}' SV_CALLERS/MERGED/ALL/delly.merged.min500bp.1suppcall.100bp_breakpoint.noTRA.vcf | sort | uniq | sed 's/<//g' | sed 's/>//g'`; do 
        if [ ! -e SV_CALLERS/MERGED/ALL/${svtype} ]; then mkdir SV_CALLERS/MERGED/ALL/${svtype}; fi;  
        # Extract SVs for class
        awk -v val=${svtype} '$5~val || $1~"#" {print}' SV_CALLERS/MERGED/ALL/delly.merged.min500bp.1suppcall.100bp_breakpoint.noTRA.vcf | \
            vcf-sort |
            sed 's/hereford\.//g' > SV_CALLERS/MERGED/ALL/${svtype}/${svtype}.vcf    
done
