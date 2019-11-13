#!/bin/bash

# Convert to vg-compliant archive
for i in {1..29}; do 
        cd CHR${i}
        hal2vg --noAncestors --hdf5InMemory --refGenome hereford CHR${i}.hal > rawCHR${i}.vg 
        cd ..; 
        echo "Done CHR ${i}" 
done

# Zip fasta folder
for i in {1..29}; do 
        cd CHR${i}
        zip -9 -r Data${i}.zip Data/ && rm -r ./Data/
        echo "Zipped chromosome $i"
        cd ..
done

