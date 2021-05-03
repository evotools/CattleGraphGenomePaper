#!/bin/bash
intervals=$1
sequences=$2
bname=`basename -s '.bed' $intervals`

# Extract upper-case sequences
bedtools getfasta -fi $sequences -bed ${intervals} -tab 2> getfasta.err | \
    python -c "import sys; [sys.stdout.write(f'{line.strip().split()[0]}\t{sum([int(c.islower()) for c in line.strip().split()[1]])}\t{len(line.strip().split()[1])}\t{line.strip().split()[1]}\n') for line in sys.stdin]" > node_analysis/regions.txt
# Add info about masked bases and sequence in the region
python 07B-combine.py ${intervals} node_analysis/regions.txt > ${bname}.masked.bed

