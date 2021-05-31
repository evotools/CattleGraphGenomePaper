#!/bin/bash
assembly=$1
reads=$2
nrun=$3

minimap2 -x map-pb -t${NSLOTS} ${assembly} ${reads} | gzip -1 > ./overlaps.${nrun}.paf.gz
