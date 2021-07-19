#!/bin/bash
intervals=$1
name=$2
outf=$3
repetitiveness=$4
#threshold=$4

module load R/3.5.3
module load roslin/samtools
module load roslin/gcc
module load anaconda
export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/minimap2
source activate bdsg
# Create folder
if [ ! -e $outf ]; then mkdir $outf; fi

# Apply filters
#awk -v var=${threshold} 'NR==1{print};NR>1 && $14<=var {print}' $intervals > ${outf}/${name}.rep${threshold}.candidate.bed
echo "" | awk 'BEGIN{OFS="\t"}; {print "SEQID", "BPI", "BPE", "NODES", "N_NODES", "STRANDS", "NODE_SEQUENCE", "N_CLOSE_TO_GAPS", "NODES_LENGTH", "REGION_SIZE", "CLASSIFICATION", "N_MASKED", "N_NT", "RATIO_MASKED", "ZSCORE", "PVAL", "SEQUENCE"}' > ${outf}/${name}.lowrep.candidate.bed
Rscript 09A-FilterRepetitive.R $intervals $repetitiveness ${outf}/${name}.lowrep.candidate.bed

# Start post processing to simplify
python -c 'import sys; [sys.stdout.write( f">{line.strip().split()[0]}_{line.strip().split()[1]}-{line.strip().split()[2]}\n{line.strip().split()[-1]}\n" ) for line in open(sys.argv[1]) if "SEQID" not in line]' ${outf}/${name}.lowrep.candidate.bed > ${outf}/${name}.lowrep.candidate.fa
samtools faidx ${outf}/${name}.lowrep.candidate.fa

source deactivate bdsg && source activate DataPy
module load igmm/apps/last

minimap2 -x asm5 --cs=long ${outf}/${name}.lowrep.candidate.fa ${outf}/${name}.lowrep.candidate.fa | \
    paftools.js view -f maf - > ${outf}/${name}.lowrep.candidate_alignments.maf
maf-convert blasttab ${outf}/${name}.lowrep.candidate_alignments.maf > ${outf}/${name}.lowrep.candidate_alignments.tmp

source deactivate DataPy && source activate bdsg && module unload igmm/apps/last
python 08B-AddScoresToBlast6.py ${outf}/${name}.lowrep.candidate_alignments.maf ${outf}/${name}.lowrep.candidate_alignments.tmp > ${outf}/${name}.lowrep.candidate_alignments.blasttab
Rscript 08C-DetectDuplicateContigs.R ${outf}/${name}.lowrep.candidate_alignments.blasttab ${outf}/${name}.lowrep.candidate.fa.fai ${outf}/${name}.lowrep.candidate.txt
python 09D-faiToBed.py ${outf}/${name}.lowrep.candidate.txt > ${outf}/${name}.lowrep.candidate.clump.bed

#python SCRIPTS/ExtractFastaSequences.py CONTIGS/contigs.nonRef.raw.fasta CONTIGS/keptContigs.noDups.txt > CONTIGS/contigs.nonRef.final.fasta