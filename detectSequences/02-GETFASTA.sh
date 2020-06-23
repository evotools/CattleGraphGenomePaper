#!/bin/bash
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N depth2seq
#$ -cwd
#$ -l h_rt=08:59:59
#$ -pe sharedmem 1
#$ -R y
#$ -l h_vmem=8.0G
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh

# Input parameters
wsize=$( python -c "import sys; print(sys.argv[-2])" $@)
minsize=$( python -c "import sys; print(sys.argv[-1])" $@)
genome=$( python -c "import sys; idx=int(sys.argv[-1]); print(sys.argv[idx])" $@ $SGE_TASK_ID)
genomes=$( python -c "import sys; print(' '.join(sys.argv[1:-2]))" $@ )
suppvec=$( python -c "import sys;tgt=sys.argv[1]; vec=['0' if i!=tgt else '1' for i in sys.argv[2:]]; print('{}#'.format(''.join(vec)))" ${genome} ${genomes} )

echo "Inputs:"
echo "All genomes: ${genomes}"
echo "Target genomes: ${genome}"
echo "Support vector: $suppvec"
echo "Distance to clump regions: ${wsize}"
echo "Minimum region size to extract: ${minsize}"

# Generate folder structure
if [ ! -e REGIONS ]; then mkdir REGIONS; fi
if [ ! -e REGIONS/${genome} ]; then mkdir REGIONS/${genome}; fi
if [ ! -e SEQUENCES/${genome} ]; then mkdir SEQUENCES/${genome}; fi

# Get the target support vector
grep $suppvec GAPLESS/${genome}/SupVect.${genome}.gapless.bed > REGIONS/${genome}/SupVect.${genome}.gapless.specific.bed
# Combine regions closer than window size
bedtools merge -d ${wsize} -i REGIONS/${genome}/SupVect.${genome}.gapless.specific.bed > REGIONS/${genome}/SupVect.${genome}.gapless.specific.concat${wsize}bp.bed;
awk -v msz=$minsize '($3-$2)>msz{print}' REGIONS/${genome}/SupVect.${genome}.gapless.specific.concat${wsize}bp.bed > REGIONS/${genome}/SupVect.${genome}.gapless.specific.concat${wsize}bp.minSize${minsize}bp.bed
# Get full fasta and fix names in file
cat /exports/cmvm/eddie/eb/groups/CTLGH_GCRF/project_andrea/BovinaeGenomes/GENOMES/${genome}.fasta | sed "s/>/>${genome}./g" > SEQUENCES/${genome}/${genome}.fasta
samtools faidx SEQUENCES/${genome}/${genome}.fasta
# Fix bed names eventually (skip if output gives error)
./SCRIPTS/FixBedNames REGIONS/${genome}/SupVect.${genome}.gapless.specific.concat${wsize}bp.minSize${minsize}bp.bed > SEQUENCES/${genome}/${genome}.bed
# Generate new fasta file
bedtools getfasta -name -fi SEQUENCES/${genome}/${genome}.fasta -fo SEQUENCES/${genome}/${genome}.specific.fasta -bed SEQUENCES/${genome}/${genome}.bed
echo "Done ${genome}"
