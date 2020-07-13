# Variant statistics
## Script to calculate metrics on the results of the different alignments

The following folder contains the script used to generate the metrics for the different VCFs generated for the 9 samples considered.

## Data  preparation
The VCFs were generated through [BAGPIPE](https://bitbucket.org/renzo_tale/bagpipe/), and then collected in folders called:
1. FB: freeayes vcf
2. GATK: GATK vcfs
3. LINEAR: vcfs generated with vg using the linear graph genome
4. CACTUS: vcfs generated with vg using the CACTUS graph genome
5. DIVERSITY: vcfs generated with vg using the DIVERSITY (CACTUS+vcfs) graph genome

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
