#!/bin/bash
#$ -cwd
#$ -N wtpoa
#$ -pe sharedmem 8
#$ -l h_rt=96:00:00
#$ -l h_vmem=32G

# Prepare dependencies
. /etc/profile.d/modules.sh
module load roslin/samtools
module load roslin/bwa/0.7.17
wtdbg2F="/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/wtdbg2"
minimap="/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/minimap2"
ucscExe='/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/UCSC'
miniasm="/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/miniasm"

# Getting arguments
while getopts ":f:a:b:o:v" opt; do
  case $opt in
    f) inputFasta=${OPTARG};;
    a) read1=${OPTARG};;
    b) read2=${OPTARG};;
    o) outname=${OPTARG};;
    v) version=${OPTARG};;
  esac
done

# polish consensus, not necessary if you want to polish the assemblies using other tools
${minimap}/minimap2 -t ${NSLOTS} -x map-pb -a ${outname}.ctg.fa $inputFasta | samtools view -Sb - >${outname}.ctg.map.bam
samtools sort -o ${outname}.ctg.map.srt -@ ${NSLOTS} ${outname}.ctg.map.bam && rm ${outname}.ctg.map.bam
samtools view ${outname}.ctg.map.srt | ${wtdbg2F}/wtpoa-cns -t ${NSLOTS} -d ${outname}.ctg.fa -i - -fo ${outname}.ctg.1st.fa
echo "Done PB polished assembly (r1)."

# polish round 2
${minimap}/minimap2 -t ${NSLOTS} -x map-pb -a ${outname}.ctg.1st.fa $inputFasta | samtools view -Sb - >${outname}.ctg.map.bam
samtools sort -o ${outname}.ctg.map.srt -@ ${NSLOTS} ${outname}.ctg.map.bam && rm ${outname}.ctg.map.bam
samtools view ${outname}.ctg.map.srt | ${wtdbg2F}/wtpoa-cns -t ${NSLOTS} -d ${outname}.ctg.1st.fa -i - -fo ${outname}.ctg.2nd.fa
echo "Done PB polished assembly (r2), and start illumina polishing."

