#!/bin/bash

export PATH=$PATH:~/Documents/Software/vcflib/bin
# Annotate each novel vcf
for i in NDama Angus Sahiwal; do 
    # Combine samples of each breed into a joint call
    bcftools merge -l ${i}.txt -O v | sed 's/hereford.//g' | vcfbreakmulti | python SVLEN.py - | bgzip -c > $PWD/${i}_final/${i}.all.vcf.gz
    tabix -p vcf $PWD/${i}_final/${i}.all.vcf.gz
    echo "File $PWD/Ndama_final/NDama.all.vcf.gz ready."
    # # Combine, annotate sizes and keep if SVLEN>500
    # bcftools merge -l ${i}.txt -O v | sed 's/hereford.//g' | vcfbreakmulti | python SVLEN.py - 500 | bgzip -c > $PWD/${i}_final/${i}.SV500.all.vcf.gz
    # echo "File $PWD/${i}_final/${i}.SV500.all.vcf.gz ready."
    if [ ! -e ${PWD}/${i}_final/VEP ]; then mkdir ${PWD}/${i}_final/VEP; fi
    mv $PWD/${i}_final/${i}.all.vcf.gz ${HOME}/vep_data/input
    docker run -t -i -v $HOME/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep vep -i /opt/vep/.vep/input/${i}.all.vcf.gz -o stdout \
        --vcf --fork 4 --species bos_taurus --variant_class --sift b --nearest symbol --distance 200 --cache --offline \
        --dir_cache /opt/vep/.vep/ \
        --dir_plugins /opt/vep/.vep/Plugins/ | bgzip -c > ${PWD}/${i}_final/${i}.all.vep.vcf.gz && \
        rm ${HOME}/vep_data/input/${i}.all.vcf.gz && \
        mv ${HOME}/vep_data/output/* $PWD/${i}_final/VEP
    gunzip -c ${PWD}/${i}_final/${i}.all.vep.vcf.gz | python SVLEN.py - 500 | bgzip -c > ${PWD}/${i}_final/${i}.all.vep.SV500.vcf.gz
    tabix -p vcf ${PWD}/${i}_final/${i}.all.vep.SV500.vcf.gz
    echo "Done VEP for ${i}"
    bedtools intersect -header -a ${PWD}/${i}_final/${i}.all.vep.SV500.vcf.gz -b INTERSECT_BN_ASM/OMregions.bed -wa | vcf-sort | vcfuniq | bgzip -c > ${PWD}/${i}_final/${i}.all.vep.SV500.Intersected.vcf.gz
    bedtools subtract -header -a ${PWD}/${i}_final/${i}.all.vep.SV500.vcf.gz -b INTERSECT_BN_ASM/OMregions.bed -wb | vcf-sort | vcfuniq | bgzip -c > ${PWD}/${i}_final/${i}.all.vep.SV500.notIntersected.vcf.gz
done

