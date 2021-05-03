#!/bin/bash
module load roslin/gcc
module load roslin/bedtools
module load roslin/augustus/3.3.3
export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/

# Inputs
bedfile=$1
genome=$2
outpath=$3
outname=$4
protdb=$5

#Ab initio prediction
if [ ! -e ${outpath} ]; then mkdir -p ${outpath}; fi
echo "Extract contig sequences..."
bedtools getfasta -fi ${genome} -bed ${bedfile} -fo ${outpath}/${outname}.fa

echo "Run augustus prediction..."
augustus --species=human ${outpath}/${outname}.fa > ${outpath}/${outname}.abInitio.gff

echo "Extract predicted sequences..."
python 11B-ExtractGenes.py -i ${outpath}/${outname}.abInitio.gff -o ${outpath}/${outname}.abInitio

# Align to protein db provided
echo "Align to database of proteins..."
diamond blastp -d ${protdb} \
    -p 4 -k 1 --more-sensitive \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore scovhsp \
    --evalue 1e-10 --id 80 -q ${outpath}/${outname}.abInitio.complete.aa > ${outpath}/${outname}.abInitio.complete.aligned.eval1e-10.id80.tab
wc -l ${outpath}/${outname}.abInitio.complete.aligned.eval1e-10.id80.tab

# Filter low-covered genes
awk '$NF > 80 {print}' ${outpath}/${outname}.abInitio.complete.aligned.eval1e-10.id80.tab > ${outpath}/${outname}.abInitio.complete.aligned.eval1e-10.id80.cov80.tab
wc -l ${outpath}/${outname}.abInitio.complete.aligned.eval1e-10.id80.cov80.tab
