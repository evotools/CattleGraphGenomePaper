#!/bin/bash
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N hal2maf
#$ -cwd
#$ -l h_rt=127:59:59
#$ -pe sharedmem 1
#$ -R y
#$ -l h_vmem=64.0G
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh
halArchive=$( python -c "import sys; print(sys.argv[-1])" $@)
genome=$( python -c "import sys; idx=int(sys.argv[-1]); print(sys.argv[idx])" $@ $SGE_TASK_ID)

if [ ! -e MAFS ]; then
    mkdir MAFS
fi
if [ ! -e MAFS ]; then
    mkdir MAFS/${genome}
fi

halStats --chromSizes ${genome} ${halArchive} > MAFS/${genome}Seqs.txt
while read line; do
    sequence=`echo $p | cut -d ' ' -f 1`
    hal2maf --refGenome ${genome} --noAncestors --inMemory BTau_align.hal MAFS/${genome}/MAF.${genome}.${sequence}.maf --refSequence ${sequence}
    echo "Done ${genome}: $sequence"
done < MAFS/${genome}Seqs.txt
