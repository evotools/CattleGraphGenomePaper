#!/bin/bash
#$ -cwd
#$ -N rmasker
#$ -pe sharedmem 24
#$ -R y
#$ -l h_rt=71:59:59
#$ -l h_vmem=16G
#$ -P roslin_ctlgh
#$ -o LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
. /etc/profile.d/modules.sh
module load igmm/apps/BEDTools
module load roslin/blast+
module load igmm/apps/BEDOPS/2.4.26
module load roslin/samtools
module load anaconda
source activate DataPy

export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/RepeatMasker/Dependencies/bins
export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/HMMER/bin
TRFPATH=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/


while getopts ":f:s:m:v" opt; do
  case $opt in
    f) fastas=${OPTARG};;
    s) species=${OPTARG};;
    m) masktype=${OPTARG};;
    v) version=${OPTARG};;
  esac
done

echo ${species}
echo ${fastas}
fasta=$( sed "$SGE_TASK_ID q;d" ${fastas} | awk '{print $1}' )
echo "Input fasta: ${fasta}"
bname=$(python -c "import sys;print sys.argv[1].strip().replace('.gz', '').replace('.fasta', '').replace('.fa', '').replace('.fna', '')" ${fasta})
echo "File base name: ${bname}"
nthreads=$(python -c "import sys; print int(sys.argv[1])/2" ${NSLOTS})
dirname=`python -c "import sys;print sys.argv[1].strip().split('/')[-1]" ${bname}`

if [ ! -e ${dirname} ]; then
        mkdir ${dirname}
fi

# Run Dustmasker and Windowmasker
dustmasker -in ${fasta} -infmt fasta -parse_seqids -outfmt interval -out ${bname}.dust.txt
awk 'BEGIN { OFS='\t'; sequence="" } { if (match($0, ">.*")) { gsub(">lcl|", "", $0); gsub(" .*", "", $0); sequence=$0 } else { gsub(" - ", "\t", $0); print sequence "\t" $1 "\t" $2+1 "\t dustmasker" } }' ${bname}.dust.txt | sed -e 's/|//g' > ${dirname}/${dirname}.dust.bed
windowmasker -in ${fasta} -infmt fasta -mk_counts -parse_seqids -sformat obinary -out ${bname}.counts
windowmasker -in ${fasta} -infmt fasta -ustat ${bname}.counts -outfmt interval -parse_seqids -out ${bname}.windowmasker.txt
awk 'BEGIN { OFS='\t'; sequence="" } { if (match($0, ">.*")) { gsub(">lcl|", "", $0); gsub(" .*", "", $0); sequence=$0 } else { gsub(" - ", "\t", $0); print sequence "\t" $1 "\t" $2+1 "\t windowsmasker"} }' ${bname}.windowmasker.txt | sed -e 's/|//g' > ${dirname}/${dirname}.windowmasker.bed

# Run RepeatMasker and combine masked regions in a bed file
/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/RepeatMasker/RepeatMasker -gff -species ${species} -dir ${dirname} -pa ${nthreads} ${fasta}
/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/RepeatMasker/util/rmOutToGFF3.pl ${dirname}/*.out > ${dirname}/${dirname}.rmask.gff
rmsk2bed < ${dirname}/*.out > ${dirname}/${dirname}.rmask.bed
cat ${dirname}/${dirname}.windowmasker.bed ${dirname}/${dirname}.dust.bed ${dirname}/${dirname}.rmask.bed | cut -f 1,2,3,4 > ${dirname}/${dirname}.all.bed.bck

samtools faidx ${dirname}/${dirname}.fasta.masked
python FixBed.py ${dirname}/${dirname}.fasta.masked.fai ${dirname}/${dirname}.all.bed.bck > ${dirname}/${dirname}.all.bed

# Perform masking
if [ $masktype == "hard" ]; then
  bedtools maskfasta -fi ${fasta} -bed ${dirname}/${dirname}.all.bed -fo ${fasta}.hardMasked.fasta
else
  bedtools maskfasta -soft -fi ${fasta} -bed ${dirname}/${dirname}.all.bed -fo ${fasta}.softMasked.fasta
fi
echo "All Done"