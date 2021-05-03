#!/bin/bash
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


while getopts ":f:s:n:o:v" opt; do
  case $opt in
    f) fasta=${OPTARG};;
    s) species=${OPTARG};;
    n) nthreads=${OPTARG};;
    o) outdir=${OPTARG};;
    v) version=${OPTARG};;
  esac
done

echo ${species}
echo ${fasta}

# Run RepeatMasker and combine masked regions in a bed file
/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/RepeatMasker/RepeatMasker -gff -species ${species} -dir ${outdir} -pa ${nthreads} ${fasta}
echo "All Done"