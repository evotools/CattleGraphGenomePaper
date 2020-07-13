#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N AB
#$ -cwd
#$ -l h_rt=4:00:00
#$ -R y
#$ -l h_vmem=1.0G
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh
#$ -q staging

. /etc/profile.d/modules.sh
module load roslin/bcftools
module load igmm/apps/tabix
module load igmm/apps/vcftools


sample=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $1}'`
vcf=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $2}'`
isectgt=`realpath $2`
tgtfld=$3
outfld=$4
currpath=`pwd`

if [ ! -e SOFT_FILT ]; then mkdir SOFT_FILT; fi
if [ ! -e SOFT_DEPTH_FILT ]; then mkdir SOFT_DEPTH_FILT; fi
if [ ! -e SOFT_DEPTH_QUAL_FILT ]; then mkdir SOFT_DEPTH_QUAL_FILT; fi

if [ ! -e $tgtfld ]; then mkdir $tgtfld; fi
if [ ! -e ${tgtfld}/$sample ]; then mkdir $tgtfld/$sample; fi
for i in SOFT_FILT SOFT_DEPTH_FILT SOFT_DEPTH_QUAL_FILT; do
    if [ ! -e ${i}/${outfld}_NOVEL ]; then mkdir ${i}/${outfld}_NOVEL; fi
    if [ ! -e ${i}/${outfld}_NOVEL/$sample ]; then mkdir ${i}/${outfld}_NOVEL/$sample; fi
    if [ ! -e ${i}/${outfld}_KNOWN ]; then mkdir ${i}/${outfld}_KNOWN; fi
    if [ ! -e ${i}/${outfld}_KNOWN/$sample ]; then mkdir ${i}/${outfld}_KNOWN/$sample; fi
    if [ ! -e ${i}/${outfld}_ALL ]; then mkdir ${i}/${outfld}_ALL; fi
    if [ ! -e ${i}/${outfld}_ALL/$sample ]; then mkdir ${i}/${outfld}_ALL/$sample; fi
done

#cd $tgtfld/$sample
bcftools isec -O z -p isec_output_$sample $vcf $isectgt && isec_output_$sample/0001.vcf.gz isec_output_$sample/0003.vcf.gz
cd ../

bcftools view -O v -f .,PASS $vcf | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT\t%DP\t%AD]\n' | \
        python ./NBySizeByQUAL.py -i - | bgzip -c > SOFT_FILT/${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.tmp.gz
    
bcftools view -O v -f .,PASS $tgtfld/$sample/isec_output_$sample/0000.vcf.gz | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT\t%DP\t%AD]\n' |\
        python ./NBySizeByQUAL.py -i - | bgzip -c > SOFT_FILT/${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.tmp.gz

bcftools view -O v -f .,PASS $tgtfld/$sample/isec_output_$sample/0002.vcf.gz | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT\t%DP\t%AD]\n' |\
        python ./NBySizeByQUAL.py -i - | bgzip -c > SOFT_FILT/${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.tmp.gz

