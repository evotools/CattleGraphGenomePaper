#!/bin/bash

CPG=$1
cd PANGENOME
if [ ! -e ./CROSSREF ]; then mkdir ./CROSSREF/; fi

for i in HF3457_Bcell_ATAC ND230_Bcell_ATAC Nelore1_Bcell_ATAC; do
        if [ ! -e ./CROSSREF/${i} ]; then
                mkdir ./CROSSREF/${i};
        fi;

        python SCRIPTS/ExtractBed.py GENRICH/${i}/*.noblank.lowrep.withQ.narrowPeak CROSSREF/${i}/TOCROSSREF.noblank.lowrep.withQ
        python SCRIPTS/SplitByGenome.py CROSSREF/${i}/TOCROSSREF.noblank.lowrep.withQ.narrowRegion.bed CROSSREF/${i}/narrow.noblank.lowrep
        echo "Validanting..."
        for w in hereford angus ankole brahman ndama; do
                upname=`python -c "import sys; print(sys.argv[1].upper())" ${w}`
                bedtools intersect -a CROSSREF/${i}/narrow.noblank.lowrep.${w}.bed -b ../${upname}/GENRICH/${i}/${i}.bed -f 0.9 -wa |\
                        sort | uniq > CROSSREF/${i}/narrow.noblank.lowrep.${w}.crossreferenced.bed
                bedtools intersect -a CROSSREF/${i}/narrow.noblank.lowrep.${w}.bed -b ../${upname}/GENRICH/${i}/${i}.bed -f 0.9 -wa -wb |\
                        sort | uniq > CROSSREF/${i}/narrow.noblank.lowrep.${w}.crossreferenced.large.bed
                echo "In CpG..."
                bedtools intersect -a CROSSREF/${i}/narrow.noblank.lowrep.${w}.crossreferenced.bed -b $CPG/CPG.${w}.bed -wa -wb > CROSSREF/${i}/narrow.noblank.lowrep.${w}.crossreferenced.inCpG.bed
        done
        cat CROSSREF/${i}/narrow.noblank.lowrep.*.crossreferenced.bed > CROSSREF/${i}/narrow.noblank.lowrep.crossreferenced.bed
        cat CROSSREF/${i}/narrow.noblank.lowrep.*.crossreferenced.large.bed > CROSSREF/${i}/narrow.noblank.lowrep.crossreferenced.large.bed
        cat CROSSREF/${i}/narrow.noblank.lowrep.*.crossreferenced.inCpG.bed > CROSSREF/${i}/narrow.noblank.lowrep.crossreferenced.inCpG.bed
        # Generate narrowPeak crossreferenced file
        echo "Prepare narrowPeak files..."
        python SCRIPTS/getNarrowPeaks.py GENRICH/${i}/*.noblank.withQ.narrowPeak CROSSREF/${i}/narrow.noblank.lowrep.crossreferenced.bed > CROSSREF/${i}/narrow.noblank.lowrep.crossreferenced.narrowPeak
        python SCRIPTS/getNarrowPeaks.py GENRICH/${i}/*.noblank.withQ.narrowPeak CROSSREF/${i}/narrow.noblank.lowrep.crossreferenced.InCpG.bed > CROSSREF/${i}/narrow.noblank.lowrep.crossreferenced.InCpG.narrowPeak
        echo "Cross-references  :      "`wc -l CROSSREF/${i}/narrow.noblank.lowrep.crossreferenced.narrowPeak`
        echo "In CpG            :      "`wc -l CROSSREF/${i}/narrow.noblank.lowrep.crossreferenced.InCpG.narrowPeak`
        # Cleaning
        echo "Cleaning..."
        for w in hereford angus ankole brahman ndama; do rm CROSSREF/${i}/*.${w}.*; done
        echo "Done $i"
        echo
done
cd ../