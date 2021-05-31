#!/bin/bash
#$ -cwd
#$ -N lrgapcl
#$ -pe sharedmem 16
#$ -l h_rt=24:00:00
#$ -l h_vmem=8G
#$ -P roslin_ctlgh

. /etc/profile.d/modules.sh
export PATH=$PATH:/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/LR_Gapcloser_v1.1


while getopts ":i:l:n" opt; do
  case $opt in
    i) inputfa=${OPTARG};;
    l) longreads=${OPTARG};;
    n) nversion=${OPTARG};;
  esac
done



/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/LR_Gapcloser_v1.1/LR_Gapcloser.sh -i ${inputfa} -l ${longreads} -s p -t ${NSLOTS}

