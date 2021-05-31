#!/bin/bash

inreads=$1
inctg=$2

refname=`basename $inctg ".fasta"`
# Prepare fasta
perl ./Solve3.3_10252018/Pipeline/10252018/fa2cmap_multi_color.pl -i ${refname}.fasta -e DLE-1 1

# Prepare the reads
filter_SNR_dynamic.pl -i ${inreads}.bnx -o ${inreads}.filter.bnx

# Run the optical mapping workflow
python ./Solve3.3_10252018/Pipeline/10252018/pipelineCL.py \
        -T 16 -j 16 -N 4 -f 0.2 -i 5 -y \
        -b ${inreads}.filter.bnx \
        -l ./output \
        -t Solve3.3_10252018/RefAligner/7915.7989rel \
        -a Solve3.3_10252018/RefAligner/7915.7989rel/optArguments_nonhaplotype_DLE1_saphyr.xml \
        -r ${refname}_DLE1_0kb_0labels.cmap \
        -C clusterConf_custom_V2.xml \
        #-B 1 -B 2 -B 3 -B 4 -B 5 -B 6 -B 7 -B 8 -B 9 -B 10 -B 11 -B 12 -B 13 -B 14 -B 15 -B 16 -B 17 -B 18 -B 19 #-B 20 -B 21 -B 22 -B 23 -B 24 #-B 25