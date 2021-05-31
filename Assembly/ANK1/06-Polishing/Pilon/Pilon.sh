#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N pilonVC
#$ -cwd
#$ -l h_rt=200:00:00
#$ -pe sharedmem 30
#$ -R y
#$ -l h_vmem=26.0G

. /etc/profile.d/modules.sh
module load java
module load roslin/bwa/0.7.17
module load roslin/samtools
module load roslin/bamtools
module load anaconda
module load igmm/apps/hdf5/1.8.16 
module load igmm/libs/glib/2.47.5
module load igmm/compilers/gcc/5.5.0

# Activate virtual environments
source activate falconpy
source /exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/FALCON/FALCON-integrate/env.sh 
source /exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/FALCON/FALCON-unzip/fc_env/bin/activate

pilon="pilon-1.23.jar"

while getopts ":a:i:p:d:o:n" opt; do
  case $opt in
    a) asm=${OPTARG};;
    i) illumina=${OPTARG};;
    p) pacbio=${OPTARG};;
    d) directory=${OPTARG};;
    o) outname=${OPTARG};;
    n) nversion=${OPTARG};;
  esac
done

if [ ${pacbio} != "" ]; then
  pbreads="--pacbio "${pacbio}
fi



java -jar -Xmx700G ${pilon} --frags ${illumina} \
                        --outdir ${directory} \
                        --output ${outname} \
                        --diploid --threads ${NSLOTS} \
                        --genome ${asm} \
                        --fix all \
                        ${pbreads}

echo "Bye."
