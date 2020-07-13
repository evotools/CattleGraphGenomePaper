#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N AB
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=1.0G
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh
#$ -q staging
. /etc/profile.d/modules.sh
module load roslin/bcftools
module load igmm/apps/tabix


sample=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $1}'`
vcf=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $2}'`
isectgt=`realpath $2`
tgtfld=$3
outfld=$4
currpath=`pwd`

if [ ! -e $tgtfld ]; then mkdir $tgtfld; fi
if [ ! -e ${tgtfld}/$sample ]; then mkdir $tgtfld/$sample; fi
if [ ! -e ${outfld}_NOVEL ]; then mkdir ${outfld}_NOVEL; fi
if [ ! -e ${outfld}_NOVEL/$sample ]; then mkdir ${outfld}_NOVEL/$sample; fi
if [ ! -e ${outfld}_KNOWN ]; then mkdir ${outfld}_KNOWN; fi
if [ ! -e ${outfld}_KNOWN/$sample ]; then mkdir ${outfld}_KNOWN/$sample; fi
if [ ! -e ${outfld}_ALL ]; then mkdir ${outfld}_ALL; fi
if [ ! -e ${outfld}_ALL/$sample ]; then mkdir ${outfld}_ALL/$sample; fi

cd $tgtfld/$sample
bcftools isec -O z -p isec_output_$sample $vcf $isectgt && isec_output_$sample/0001.vcf.gz isec_output_$sample/0003.vcf.gz
cd ../../
bcftools view -f .,PASS -i 'GT="het"' $vcf | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' | \
    python ${currpath}/ABbyIndelSizeV1.2.py -i - | bgzip -c > ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.txt.gz

bcftools view -f .,PASS -i 'GT="het"' $tgtfld/$sample/isec_output_$sample/0000.vcf.gz | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' | \
    python ${currpath}/ABbyIndelSizeV1.2.py -i - | bgzip -c > ${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.txt.gz

bcftools view -f .,PASS -i 'GT="het"' $tgtfld/$sample/isec_output_$sample/0002.vcf.gz | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' | \
    python ${currpath}/ABbyIndelSizeV1.2.py -i - | bgzip -c > ${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.txt.gz


