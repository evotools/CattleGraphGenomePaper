#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N AB
#$ -cwd
#$ -l h_rt=6:00:00
#$ -l h_vmem=1.0G
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh
#$ -q staging
. /etc/profile.d/modules.sh
module load roslin/bcftools
module load igmm/apps/tabix
module load igmm/apps/vcftools
module load igmm/apps/plink
export PATH=$PATH:


vcffile=$1
isectgt=`realpath $2`
tgtfld=$3
outfld=$4
currpath=`pwd`
ALG_NAME=`python -c "import sys; print(sys.argv[1].split('/')[-1])" $tgtfld`

if [ ! -e $outfld ]; then mkdir $outfld; fi
if [ ! -e $outfld/${ALG_NAME} ]; then mkdir $outfld/${ALG_NAME}; fi

# Combine all vcfs
echo "Merging"
awk '{print $2}' $vcffile > ${vcffile}.list
bcftools merge -f .,PASS -l $vcffile.list -O v | \
    vcftools --vcf - --non-ref-ac 1 --recode --recode-INFO-all --stdout | \
    bcftools view -O b -o $outfld/${ALG_NAME}/${ALG_NAME}.bcf 
echo "Indexing"
bcftools index $outfld/${ALG_NAME}/${ALG_NAME}.bcf 
# Intersect with known sites
echo "Intersecting"
bcftools isec -O b -p $outfld/${ALG_NAME}/isec_output $outfld/${ALG_NAME}/${ALG_NAME}.bcf $isectgt \
    && rm $outfld/${ALG_NAME}/isec_output/0001.bcf $outfld/${ALG_NAME}/isec_output/0003.bcf

# Extract missingness
echo "Computing missingness"
plink --cow --allow-extra-chr --bcf $outfld/${ALG_NAME}/${ALG_NAME}.bcf --missing --out $outfld/${ALG_NAME}/${ALG_NAME}.missing
plink --cow --allow-extra-chr --bcf $outfld/${ALG_NAME}/isec_output/0000.bcf --missing --out $outfld/${ALG_NAME}/${ALG_NAME}.novel.missing
plink --cow --allow-extra-chr --bcf $outfld/${ALG_NAME}/isec_output/0002.bcf --missing --out $outfld/${ALG_NAME}/${ALG_NAME}.known.missing