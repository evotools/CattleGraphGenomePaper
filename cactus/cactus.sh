#!/bin/bash

# Create data folder and structuree
for i in {1..29}; do 
        mkdir ./CHR${i}/
	mkdir ./CHR${i}/Data
        for j in hereford ndama angus brahman ankole; do 
                samtools faidx ../Genomes/${j}.fasta $i > CHR${i}/Data/${j}${i}.fasta
        done
        sed "s/CHR/$i/g" conffile.txt > CHR${i}/conffile.txt
        echo "Done ${i}"
done

# Run cactus
for i in {1..29}; do 
        cd CHR${i}
        /usr/bin/time -v cactus jobStore conffile.txt CHR${i}.hal --maxCores 16 --maxMemory 96G
        cd ..
done

