#!/bin/bash
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N getNotRef
#$ -cwd
#$ -l h_rt=08:59:59
#$ -pe sharedmem 1
#$ -R y
#$ -l h_vmem=8.0G
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh
. /etc/profile.d/modules.sh
module load roslin/gcc
module load igmm/libs/boost
module load roslin/bedtools
module load igmm/apps/last
module load R/3.5.3

export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/silix-1.2.9/src
export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/minimap2/


# Input parameters
wsize=$( python -c "import sys; print(sys.argv[-2])" $@)
minsize=$( python -c "import sys; print(sys.argv[-1])" $@)
genome=$( python -c "import sys; idx=int(sys.argv[-1]); print(sys.argv[idx])" $@ $SGE_TASK_ID)
genomes=$( python -c "import sys; idx=int(sys.argv[-1]); [sys.stdout.write('{} '.format(i)) for n,i in enumerate(sys.argv[1:-1]) if n+1 != idx ]" $@ $SGE_TASK_ID)
#genomes=$( python -c "import sys; print(' '.join(sys.argv[1:-2]))" $@ )

echo "Inputs:"
echo "Target genomes: ${genomes}"
echo "Removed genome: ${genome}"
echo "Distance to clump regions: ${wsize}"
echo "Minimum region size to extract: ${minsize}"

# Remove hereford 
python SCRIPTS/RemoveOneGenome.py $( for genome in ${genomes}; do echo GAPLESS/${genome}/SupVect.${genome}.gapless.bed; done ) $SGE_TASK_ID


# Generate folder structure
for genome in ${genomes}; do
    if [ ! -e NONREF ]; then mkdir NONREF; fi
    if [ ! -e NONREF/${genome} ]; then mkdir NONREF/${genome}; fi
    mv GAPLESS/${genome}/SupVect.${genome}.gapless.nonRef.bed NONREF/${genome}

    # Combine nearby regions
    bedtools merge -d ${wsize} -i NONREF/${genome}/SupVect.${genome}.gapless.nonRef.bed > NONREF/${genome}/SupVect.${genome}.gapless.nonRef.concat${wsize}bp.bed;

    # Extract by minimum size provided
    awk -v msz=$minsize '($3-$2)>msz{print}' NONREF/${genome}/SupVect.${genome}.gapless.nonRef.concat${wsize}bp.bed > NONREF/${genome}/SupVect.${genome}.gapless.nonRef.concat${wsize}bp.minSize${minsize}bp.bed

    ## Get full fasta and fix names in file
    # Fix bed names eventually (skip if output gives error)
    ./SCRIPTS/FixBedNames NONREF/${genome}/SupVect.${genome}.gapless.nonRef.concat${wsize}bp.minSize${minsize}bp.bed > NONREF/${genome}/${genome}.contigs.bed
    # Generate new fasta file
    bedtools getfasta -name -fi SEQUENCES/${genome}/${genome}.fasta -fo NONREF/${genome}/${genome}.contigs.nonRef.fasta -bed NONREF/${genome}/${genome}.contigs.bed
    echo "Done ${genome}"
done

cat NONREF/*/*.contigs.nonRef.fasta > CONTIGS/contigs.nonRef.raw.fasta
samtools faidx CONTIGS/contigs.nonRef.raw.fasta

# Self script
minimap2 -x asm5 --cs=long CONTIGS/contigs.nonRef.raw.fasta CONTIGS/contigs.nonRef.raw.fasta | \
    paftools.js view -f maf - > CONTIGS/contigs.maf
maf-convert blasttab CONTIGS/contigs.maf > CONTIGS/contigs.tmp
python SCRIPTS/AddScoresToBlast6.py CONTIGS/contigs.maf CONTIGS/contigs.tmp > CONTIGS/contigs.blasttab
SCRIPTS/DetectDuplicateContigs.R CONTIGS/contigs.blasttab CONTIGS/contigs.nonRef.raw.fasta.fai CONTIGS/keptContigs.noDups.txt
python SCRIPTS/ExtractFastaSequences.py CONTIGS/contigs.nonRef.raw.fasta CONTIGS/keptContigs.noDups.txt > CONTIGS/contigs.nonRef.final.fasta