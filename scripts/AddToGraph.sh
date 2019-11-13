#!/bin/bash
#
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N graphAdd
#$ -cwd
#$ -l h_rt=47:00:00
#$ -pe sharedmem 1
#$ -R y
#$ -l h_vmem=48.0G
#$ -hold_jid gc_60419
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh
. /etc/profile.d/modules.sh
module load java/jdk/1.8.0
module load roslin/samtools/1.9
module load igmm/compilers/gcc
module load igmm/apps/tabix
export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/vg


# 
# Start program
# 
while getopts ":g:v:s:n" opt; do
  case $opt in
    g) gv=${OPTARG};;
    v) vc=${OPTARG};;
    s) sn=${OPTARG};;
    n) nversion=${OPTARG};;
  esac
done




graph=`head -n $SGE_TASK_ID $gv | tail -1`
if [ ! -z $vc ]; then
	vcf=`head -n $SGE_TASK_ID $vc | tail -1`
	spp=$sn
	if [ ! -e $vcf".tbi" ]; then tabix -p vcf $vcf; fi
fi
gName=`basename -s ".vg" $graph`
export TMPDIR=`pwd`

echo $vcf
echo $spp
echo $graph

if [ ! -z $vcf ] && [ ! -z $spp ] ; then 
	echo $spp
	echo "${vg}/vg add -n ${SGE_TASK_ID}=${spp}.${SGE_TASK_ID} -v $vcf $graph > GRAPH/${gName}.tmp.vg"
	vg mod -X 32 $graph | vg ids -s - | vg add -n ${SGE_TASK_ID}=${spp}.${SGE_TASK_ID} -v $vcf > GRAPH/${gName}.tmp.vg
	vg ids -s GRAPH/${gName}.tmp.vg > GRAPH/${gName}.final.vg && rm GRAPH/${gName}.tmp.vg
else
	vg mod -X 32 $graph | vg ids -s - > GRAPH/${gName}.final.vg
fi
