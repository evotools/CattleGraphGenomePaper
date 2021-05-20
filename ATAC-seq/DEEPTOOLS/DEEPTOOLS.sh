#!/bin/bash
#$ -N deeptools
#$ -R y
#$ -pe sharedmem 4
#$ -l h_vmem=8G
#$ -l h_rt=23:59:59
#$ -P roslin_ctlgh
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.err
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.out
#$ -cwd

. /etc/profile.d/modules.sh
module load anaconda
source activate deeptools

echo "Compute general matrix and plot"
computeMatrix scale-regions -S DEEPTOOLS/MULTI/HF3457_Bcell_ATAC/HF3457_Bcell_ATAC.bw DEEPTOOLS/MULTI/ND230_Bcell_ATAC/ND230_Bcell_ATAC.bw DEEPTOOLS/MULTI/Nelore2_Bcell_ATAC/Nelore2_Bcell_ATAC.bw -R FEATURES/Bos_taurus.ARS-UCD1.2.103.chr.gtf --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 -o ./DEEPTOOLS/MULTI/matrix.genes.multi.mat.gz --missingDataAsZero -p 4 --skipZeros --samplesLabel "Holstein 1" "N'Dama 2" "Nelore 1"
plotHeatmap -m ./DEEPTOOLS/MULTI/matrix.genes.multi.mat.gz -out ./DEEPTOOLS/MULTI/multi_GenesHeatmap.png

echo "Compute novel contigs matrix and plot"
computeMatrix scale-regions -S DEEPTOOLS/MULTI/HF3457_Bcell_ATAC/HF3457_Bcell_ATAC.bw DEEPTOOLS/MULTI/ND230_Bcell_ATAC/ND230_Bcell_ATAC.bw DEEPTOOLS/MULTI/Nelore2_Bcell_ATAC/Nelore2_Bcell_ATAC.bw -R FEATURES/predicted.abInitio.bed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 -o ./DEEPTOOLS/MULTI/matrix.augustus.multi.mat.gz --missingDataAsZero -p 4 --skipZeros --samplesLabel "Holstein 1" "N'Dama 2" "Nelore 1"
plotHeatmap -m ./DEEPTOOLS/MULTI/matrix.augustus.multi.mat.gz -out ./DEEPTOOLS/MULTI/multi_augustusGenesHeatmap.png

# repeat reference point
computeMatrix reference-point -S DEEPTOOLS/MULTI/HF3457_Bcell_ATAC/HF3457_Bcell_ATAC.bw DEEPTOOLS/MULTI/ND230_Bcell_ATAC/ND230_Bcell_ATAC.bw DEEPTOOLS/MULTI/Nelore2_Bcell_ATAC/Nelore2_Bcell_ATAC.bw -R FEATURES/Bos_taurus.ARS-UCD1.2.103.chr.gtf -a 3000 -b 3000 -o ./DEEPTOOLS/MULTI/matrix.genes.multi.refpoint.mat.gz --missingDataAsZero -p $NSLOTS --skipZeros --samplesLabel "Holstein 1" "N'Dama 2" "Nelore 1"
plotHeatmap -m ./DEEPTOOLS/MULTI/matrix.genes.multi.refpoint.mat.gz -out ./DEEPTOOLS/MULTI/multi_GenesHeatmap_refpoint.png

computeMatrix reference-point -S DEEPTOOLS/MULTI/HF3457_Bcell_ATAC/HF3457_Bcell_ATAC.bw DEEPTOOLS/MULTI/ND230_Bcell_ATAC/ND230_Bcell_ATAC.bw DEEPTOOLS/MULTI/Nelore2_Bcell_ATAC/Nelore2_Bcell_ATAC.bw -R FEATURES/predicted.abInitio.bed -a 3000 -b 3000 -o ./DEEPTOOLS/MULTI/matrix.augustus.multi.refpoint.mat.gz --missingDataAsZero -p $NSLOTS --skipZeros --samplesLabel "Holstein 1" "N'Dama 2" "Nelore 1"
plotHeatmap -m ./DEEPTOOLS/MULTI/matrix.augustus.multi.refpoint.mat.gz -out ./DEEPTOOLS/MULTI/multi_augustusGenesHeatmap_refpoint.png
