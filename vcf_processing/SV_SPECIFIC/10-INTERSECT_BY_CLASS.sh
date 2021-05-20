#!/bin/bash

echo "Intersect by class of SV"

if [ ! -e INTERSECTIONS ]; then mkdir INTERSECTIONS; fi
for svtype in DEL DUP INS INV COMPLEX; do
    if [ ! -e INTERSECTIONS/${svtype} ]; then mkdir INTERSECTIONS/${svtype}; fi

    # Intersect OM and VG if both have class of SV
    if [ -e JOINT/CLASSES/${svtype} ] && [ -e BNOM/BED_V2/${svtype} ]; then 
            if [ ! -e INTERSECTIONS/${svtype}/VG5P ]; then mkdir INTERSECTIONS/${svtype}/VG5P; fi
            bedtools intersect -a JOINT/CLASSES/${svtype}/${svtype}.bed -b BNOM/BED_V2/${svtype}/${svtype}.sort.bed -wa -u |\
                bedtools sort -i - > INTERSECTIONS/${svtype}/VG5P/VG5P.${svtype}.isecOM.bed
    fi 

    if [ -e SV_CALLERS/MERGED/ALL/${svtype} ] && [ -e BNOM/BED_V2/${svtype} ]; then 
            if [ ! -e INTERSECTIONS/${svtype}/DELLY ]; then mkdir INTERSECTIONS/${svtype}/DELLY; fi
            bedtools intersect -a SV_CALLERS/MERGED/ALL/${svtype}/${svtype}.vcf -b BNOM/BED_V2/${svtype}/${svtype}.sort.bed -header -wa -u |\
                vcf-sort > INTERSECTIONS/${svtype}/DELLY/DELLY.${svtype}.isecOM.vcf
    fi 
done

cat INTERSECTIONS/*/VG5P/*.bed | bedtools sort -i - > INTERSECTIONS/VG5P_isecOM.bed
vcf-concat INTERSECTIONS/*/DELLY/*.vcf | vcf-sort | vcfcreatemulti > INTERSECTIONS/DELLY_isecOM.vcf