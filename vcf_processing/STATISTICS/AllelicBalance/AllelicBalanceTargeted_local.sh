#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
. /etc/profile.d/modules.sh
module load roslin/bcftools
module load igmm/apps/tabix

FILELIST=$1
isectgt=`realpath $2`
tgtfld=$3
outfld=$4
currpath=`pwd`

while read p; do
	sample=`echo $p | awk '{print $1}'`
	vcf=`echo $p | awk '{print $2}'`
	if [ ! -e $tgtfld ]; then mkdir $tgtfld; fi
	if [ ! -e ${tgtfld}/$sample ]; then mkdir $tgtfld/$sample; fi
	cd $tgtfld/$sample
	bcftools isec -O z -p isec_output_$sample $vcf $isectgt && rm isec_output_$sample/0001.vcf.gz isec_output_$sample/0003.vcf.gz

	bcftools view -f .,PASS -i 'GT="het"' $vcf | \
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' | \
		python ${currpath}/ABbyIndelSize.py - | bgzip -c > ${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.txt.gz

	bcftools view -f .,PASS -i 'GT="het"' isec_output_$sample/000.vcf.gz | \
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' | \
		python ${currpath}/ABbyIndelSize.py - | bgzip -c > ${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.txt.gz

	bcftools view -f .,PASS -i 'GT="het"' isec_output_$sample/0002.vcf.gz | \
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' | \
		python ${currpath}/ABbyIndelSize.py - | bgzip -c > ${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.txt.gz

	cd ../../
	echo $sample

done < $FILELIST


