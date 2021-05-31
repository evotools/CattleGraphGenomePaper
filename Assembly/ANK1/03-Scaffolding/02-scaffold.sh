#!/bin/bash
module load R
module load roslin/python/2.7.13

inreads=$1
inctg=$2

# Run hybrid scaffolding
perl ./Solve3.3_10252018/HybridScaffold/12162019//hybridScaffold.pl \
        -n ${refname}.fa \
        -b ./output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap \
        -c ./hybridScaffold_DLE1_config.xml \
        -r ./Solve3.3_10252018/RefAligner/10330.10436rel/RefAligner \
        -o hybridScfld_nonHap -f -B 2 -N 2 -y \
        -m ${inreads}.filter.bnx \
        -e ./output/contigs/auto_noise/autoNoise1.errbin \
        -q optArguments_nonhaplotype_DLE1_saphyr.xml \
        -p ./Solve3.3_10252018/Pipeline/12162019/ \
        -z BionanoHybridScaffoldingNonHap.tar.gz
