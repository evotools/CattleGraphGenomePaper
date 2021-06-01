#!/bin/bash
#$ -cwd
#$ -N ragout2
#$ -pe sharedmem 4
#$ -l h_rt=24:00:00
#$ -l h_vmem=8G
#$ -P roslin_ctlgh

. /etc/profile.d/modules.sh
module load anaconda
module load roslin/ragout
source activate DataPy

export CLAPACKPATH=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/clapack
export ENABLE_PHYLOP=1
export PATH=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/hdf5/bin:${PATH}
export h5prefix=-prefix=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/hdf5
export PATH=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/hal/bin:${PATH}
export PYTHONPATH=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/:${PYTHONPATH}

ragout $1 --refine --outdir ./$2 -t $NSLOTS -s maf --overwrite