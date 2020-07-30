# Get SV for specific individuals
## Script to extract SV specific for a single breed.

The following folder contains the script used to extract the large SV specific for each breed, excluding the ones identified in the other populations considered.
Briefly, the script extract from each sample the varaints not found in any individuals of the other populations using ```bcftools isec```. The resulting variants are then compared with variants identified by the optical mapping for two samples, and converted to bed file with large interval.

## Optical mapping data preparation
After processing the bnx files through Bionano Solve pipeline, we convert the final smap to vcf using the ```smap2vcf_v2.py``` script. Then, we removed the translocations and converted variants to  using the command:
```
bcftools query -f'%CHROM\t%POS\t%END\t\t%ID#%SVTYPE#LEN=%SVLEN#CIEND=%CIEND#CIPOS=%CIPOS\n' sample.vcf.gz | \ 
    grep -v "#TRA#" | \ 
    sort -k1,1n -k2,2n -k3,3n > ../BED/sample.SV.noTranslocations.bed 
```
The result for each sample are then concatenated and sorted:
```
mkdir INTERSECT_BN_ASM
cat sample{1..2}.SV.noTranslocations.bed | bedtools sort -i - > INTERSECT_BN_ASM/OMregions.bed
```


# Processing the VCF files


Withing each folder there should be a subfolder with the sample name containing a vcf for each sample. Files will then be listed in different text lists that will be used as input for the different analyses.
In addition, a vcf with the variants included in the cactus graph genome has been generated through deconstructing each vg graph:
```
vg snarls -t 1 GRAPH_CACTUS/CHRN.final.vg > GRAPH_CACTUS/CHRN.snarls
vg deconstruct -e -P hereford -A angus,ndama,ankole,brahman -r GRAPH_CACTUS/CHRN.snarls -t 1 GRAPH_CACTUS/CHRN.final.vg | bgzip -c > DECONSTRUCTED/CHRN.vcf.gz
```
Where N is the number of the chromosome to be processed. Then, vcfs have been concatenated using vcf-concat:
```
vcf-concat CHR{1..29}.vcf.gz | bgzip -c > DECONSTRUCTED/CACTUSVAR.vcf.gz
```

## Run the single pipelines
Within each folder there is a series of multiqsub.sh scripts, that will submit the analyses for each algorith separately. The analyses will briefly intersect each file with the variants presents in the cactus graph, and classify them as 
1. Novel: variant not present in the cactus graph
2. Known: variants in the cactus graph
3. All: the complete set of variants 

After extracting the parameter of interest (allelic balance, TsTv, Size, etc), it will proceed summarising the results by variant size using custom python scripts.
Finally, the plots can be generated through the scripts in the folder ```PLOT_R```.
