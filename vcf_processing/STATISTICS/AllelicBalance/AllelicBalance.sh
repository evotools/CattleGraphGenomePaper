#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N AB
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=2.0G
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh
. /etc/profile.d/modules.sh
module load roslin/bcftools
module load igmm/apps/tabix


sample=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $1}'`
vcf=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $2}'`
isectgt=`realpath $2`
tgtfld=$3
currpath=`pwd`
if [ ! -e $tgtfld ]; then mkdir $tgtfld; fi
if [ ! -e ${tgtfld}/$sample ]; then mkdir $tgtfld/$sample; fi

cd $tgtfld/$sample
#bcftools isec -O z -p isec_output_$sample $vcf $isectgt 
#bcftools view -i 'GT="het"' isec_output_$sample/0002.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' | python ${currpath}/ABbyIndelSize.py - > ${sample}_AllelicBalanceKnownSites.txt

#bcftools view -i 'GT="het" & INFO/DP>=15 & QUAL>100' $vcf | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' | python ${currpath}/ABbyIndelSize.py - > ${sample}_AllelicBalanceAllSites_filt.txt
bcftools view -i 'GT="het" & INFO/DP>=15 & QUAL>100' isec_output_$sample/0000.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' | python ${currpath}/ABbyIndelSize.py - > ${sample}_AllelicBalanceNovelSites_filt.txt
bcftools view -i 'GT="het" & INFO/DP>=15 & QUAL>100' isec_output_$sample/0002.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' | python ${currpath}/ABbyIndelSize.py - > ${sample}_AllelicBalanceKnownSites_filt.txt

#bcftools view -i 'GT="het"' $vcf | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' | python ${currpath}/ABbyIndelSize.py - > ${sample}_AllelicBalanceAllSites.txt
#python ${currpath}/SummariseAB.py ${sample}_AllelicBalanceNovelSites.txt
echo "Done $sample"

