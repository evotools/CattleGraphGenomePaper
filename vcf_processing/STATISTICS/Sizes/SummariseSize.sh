bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N ABsum
#$ -cwd
#$ -l h_rt=00:20:00
#$ -R y
#$ -l h_vmem=8.0G
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh

. /etc/profile.d/modules.sh
# module load roslin/bcftools
# module load igmm/apps/tabix
module load anaconda
module load roslin/gcc
source activate DataPy
module load igmm/apps/tabix

sample=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $1}'`
vcf=`head -$SGE_TASK_ID $1 | tail -1 | awk '{print $2}'`
isectgt=`realpath $2`
tgtfld=$3
outfld=$4
currpath=`pwd`

python ${currpath}/NormaliseQUAL.py -i SOFT_FILT/${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.tmp.gz | python SummariseSize.py > SOFT_FILT/${outfld}_ALL/$sample/${sample}.csv && rm SOFT_FILT/${outfld}_ALL/$sample/${sample}_AllelicBalanceAllSites.tmp.gz
python ${currpath}/NormaliseQUAL.py -i SOFT_FILT/${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.tmp.gz | python SummariseSize.py > SOFT_FILT/${outfld}_NOVEL/$sample/${sample}.csv && rm SOFT_FILT/${outfld}_NOVEL/$sample/${sample}_AllelicBalanceNovelSites.tmp.gz
python ${currpath}/NormaliseQUAL.py -i SOFT_FILT/${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.tmp.gz | python SummariseSize.py > SOFT_FILT/${outfld}_KNOWN/$sample/${sample}.csv && tm SOFT_FILT/${outfld}_KNOWN/$sample/${sample}_AllelicBalanceKnownSites.tmp.gz

