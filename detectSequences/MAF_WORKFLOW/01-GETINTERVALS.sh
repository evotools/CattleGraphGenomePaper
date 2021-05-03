#!/bin/bash
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N CalcCov
#$ -cwd
#$ -l h_rt=47:59:59
#$ -pe sharedmem 4
#$ -R y
#$ -l h_vmem=128.0G
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh
export TMPDIR=`pwd`
. /etc/profile.d/modules.sh
module load anaconda
module load roslin/bedtools
source activate DataPy3
scriptPaths=$PWD/SCRIPTS/

genomes=`python -c "import sys; print(','.join(sys.argv[1:]))" $@`
fasta=`python -c "import sys; print(sys.argv[int(sys.argv[-1])])" $@ $SGE_TASK_ID`
echo "Groups: $genomes"
echo "Processing $fasta"

if [ ! -e SUPPORTBED ]; then mkdir SUPPORTBED; fi 
if [ ! -e SUPPORTBED/${fasta} ]; then mkdir SUPPORTBED/${fasta}; fi 
if [ ! -e SUPPORTBED/${fasta}/TMP ]; then mkdir SUPPORTBED/${fasta}/TMP; fi 
if [ ! -e GAPLESS ]; then mkdir GAPLESS; fi  
if [ ! -e GAPLESS/${fasta} ]; then mkdir GAPLESS/${fasta}; fi  
if [ ! -e SUMMARIES ]; then mkdir SUMMARIES; fi 
if [ ! -e TMP ]; then mkdir TMP; fi

# Process each chromosome separately
while read p; do 
    seqName=`echo $p | awk '{print $1}'`
    python $scriptPaths/MAFtoCoverage.py -m MAFS/${fasta}/MAF.${fasta}.${seqName}.maf -t $NSLOTS -o stdout --highmem | \
        awk -v var=${fasta}.${seqName} '$1==var {print}' | \
        python $scriptPaths/CoverageToSupport.py - ${genomes} stdout | \
        bedtools sort -i - | \
        python $scriptPaths/BigSupportToSmallSupport.py - stdout | \
        bedtools sort -i - > SUPPORTBED/$fasta/TMP/SupVect.${fasta}.${seqName}.bed
    echo "Done $bname"
done < MAFS/${fasta}Seqs.txt

# Combine the single bed files into a large final bed file
for i in `ls ./MAFS/${fasta} | grep '.maf'`; do 
    seqName=`echo $p | awk '{print $1}'`
    cat SUPPORTBED/$fasta/TMP/SupVect.${fasta}.${seqName}.bed #&& rm SUPPORTBED/$fasta/TMP/SupVect.${fasta}.${seqName}.bed
done < MAFS/${fasta}Seqs.txt > SUPPORTBED/$fasta/SupVect.${fasta}.bed

# Extract regions of interest and create an upset plot
grep -v -e "#0#1#" -e "#1#1#" SUPPORTBED/$fasta/SupVect.${i}.bed > GAPLESS/${fasta}/SupVect.${fasta}.gapless.bed
python $scriptPaths/UpsetPy.py -i GAPLESS/${fasta}/SupVect.${fasta}.gapless.bed -t bed -o SUMMARIES/SummaryNoGaps_${fasta} -g $genomes
python $scriptPaths/UpsetPy.py -n -i GAPLESS/${fasta}/SupVect.${fasta}.gapless.bed -t bed -o SUMMARIES/SummaryNoGapsNoReps_${fasta} -g $genomes
echo "Done ${fasta}"



