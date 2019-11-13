#!/bin/bash

# Generate folder
if [ ! -e GRAPH ]; then
    mkdir GRAPH
fi

# Create initial chromosome-level VG 
qsub -t 1-29 -tc 29 AddToGraph.sh -g ../LISTS/rawvgList.txt -v ../LISTS/vcflist.txt -s hereford

# Create common ID space
qsub GenerateIds.sh

# Genreate XG index
qsub GenerateXg.sh

# Generate mapping of the nodes
qsub GenerateMap.sh

# Prune graph
qsub PruneGraph.sh -g rawvgList.txt

# Generate GCSA index
qsub GenerateGCSA.sh

# Combine the single chromosome vg into one large vg
cd GRAPH
cat `for i in {1..29}; do echo CHR${i}.vg; done` > all.vg
cd ../