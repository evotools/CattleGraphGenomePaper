#!/bin/bash

for i in HF3457_Bcell_ATAC ND230_Bcell_ATAC Nelore1_Bcell_ATAC; do 
    bedtools sort -i PANGENOME_NONHEREFORD_20M/GENRICH/${i}/${i}.noblank.narrowPeak | cut -f 1-3,5 > PANGENOME_NONHEREFORD_20M/GENRICH/${i}/${i}.noblank.sort.bed
    bedGraphToBigWig PANGENOME_NONHEREFORD_20M/GENRICH/${i}/${i}.noblank.sort.bed ./PANGENOME_NONHEREFORD_20M/pangenome.genome PANGENOME_NONHEREFORD_20M/GENRICH/${i}/${i}.noblank.sort.bw
done

# Make matrix and plot for the genes, skipping the missing ones
computeMatrix scale-regions \
    -S PANGENOME_NONHEREFORD_20M/GENRICH/HF3457_Bcell_ATAC/HF3457_Bcell_ATAC.noblank.sort.bw \
        PANGENOME_NONHEREFORD_20M/GENRICH/ND230_Bcell_ATAC/ND230_Bcell_ATAC.noblank.sort.bw \
        PANGENOME_NONHEREFORD_20M/GENRICH/Nelore2_Bcell_ATAC/Nelore2_Bcell_ATAC.noblank.sort.bw \
    -R FEATURES/Bos_taurus.ARS-UCD1.2.103.chr.gtf \
    --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
    --skipZeros --missingDataAsZero \ 
    -p 4 \
    -o matrix.all_genes.mat.gz
plotHeatmap -m matrix.all_genes.mat.gz -out allGenesHeatmap.png

# Make matrix and plot for the augustus-predicted genes, skipping the missing ones
computeMatrix scale-regions \
    -S PANGENOME_NONHEREFORD_20M/GENRICH/HF3457_Bcell_ATAC/HF3457_Bcell_ATAC.noblank.sort.bw \
        PANGENOME_NONHEREFORD_20M/GENRICH/ND230_Bcell_ATAC/ND230_Bcell_ATAC.noblank.sort.bw \
        PANGENOME_NONHEREFORD_20M/GENRICH/Nelore2_Bcell_ATAC/Nelore2_Bcell_ATAC.noblank.sort.bw \
    -R FEATURES/predicted.abInitio.bed \
    --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
    --skipZeros --missingDataAsZero \ 
    -p 4 \
    -o matrix.augustus_genes.mat.gz
plotHeatmap -m matrix.augustus_genes.mat.gz -out predictedGenesHeatmap.png