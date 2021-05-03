#!/bin/bash
# Add vcflib
export PATH=$PATH:~/Documents/Software/vcflib/bin

# Extract vcf with abs(SVLEN) > 500 bp
bcftools merge -m all -l samples.txt -O v | sed 's/hereford.//g' | vcfbreakmulti | python SVLEN.py - 500 | bgzip -c > JOINT/JOINED.SV500.vcf.gz

# Get all indels
bcftools merge -m all -l samples.txt -O v | sed 's/hereford.//g' | vcfbreakmulti | bcftools view -v indels -O z -o JOINT/JOINED.indels.vcf.gz

# Run VEP with Docker
if [ ! -e ${HOME}/vep_data/input ]; then mkdir -p ${HOME}/vep_data/input; fi
mv JOINT/JOINED.indels.vcf.gz ${HOME}/vep_data/input
docker run -t -i -v $HOME/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep vep -i /opt/vep/.vep/input/JOINED.indels.vcf.gz -o stdout \
    --vcf --fork 4 --species bos_taurus --variant_class --sift b --nearest symbol --distance 200 --cache --offline \
    --dir_cache /opt/vep/.vep/ \
    --dir_plugins /opt/vep/.vep/Plugins/ | \
        bgzip -c > ${PWD}/JOINT/JOINED.vep.indels.vcf.gz && \
        rm ${HOME}/vep_data/input/JOINED.indels.vcf.gz && \
        mv ${HOME}/vep_data/output/* $PWD/JOINT/VEP

# Extract inframe and frameshift variants
gunzip -c JOINT/JOINED.vep.indels.vcf.gz| python SVLEN.py - | awk '$1~"#" || $0~"frameshift" || $0~"inframe" {print}' | bgzip -c > INDELS_SIZE_INFRAME/JOINT.frame.vep.indels.vcf.gz
