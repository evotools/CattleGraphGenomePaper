#!/bin/bash
# Script to process peaks in each file

for w in HEREFORD ANGUS NDAMA ANKOLE BRAHMAN PANGENOME; do
    for i in HF3457_Bcell_ATAC ND230_Bcell_ATAC Nelore1_Bcell_ATAC; do 
        init=`python -c "import sys; print(sum([1 for i in open(sys.argv[1]) ] ))" $PWD/${w}/GENRICH/${i}/${i}.narrowpeak`
        noblank=`python -c "import sys; print(sum([1 for i in open(sys.argv[1]) ] ))" $PWD/${w}/GENRICH/${i}/${i}.noblank.narrowpeak`
        lowrep=`python -c "import sys; print(sum([1 for i in open(sys.argv[1]) ] ))" $PWD/${w}/GENRICH/${i}/${i}.noblank.lowrep.narrowpeak`
        if [ $w == "PANGENOME" ]; then
            crossref=`python -c "import sys; print(sum([1 for i in open(sys.argv[1]) ] ))" ${w}/CROSSREF/${i}/narrow.noblank.lowrep.crossreferenced.narrowPeak`
            crossCPG=`python -c "import sys; print(sum([1 for i in open(sys.argv[1]) ] ))" ${w}/CROSSREF/${i}/narrow.noblank.lowrep.crossreferenced.InCpG.narrowPeak`
        else
            crossref=0
            crossCPG=0
        fi
    echo $w $i $init $noblank $lowrep $crossref $crossCPG
    done
done