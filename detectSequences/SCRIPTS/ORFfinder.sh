#!/bin/bash

export PATH=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea:$PATH

python ./SCRIPTS/rename.py ${1} ${2}/renamed.fasta
ORFfinder -in ${2}/renamed.fasta -s 0 -n true > ${2}/orffinder.aa