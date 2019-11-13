# GraphGenomeScripts
Scripts used for the graph genome indexing

## Running CACTUS
CACTUS has been run on single chromosomes on a [Eleanor cloud](https://www.ed.ac.uk/information-services/computing/computing-infrastructure/cloud-computing-service) instance with
16 cores, 96Gb of RAM and 320 Gb of disk space. HAL2VG has been on each chromosome separately on the same virtual machine. 
The script used to perform the alignments, as well as the general configuration file for CACTUS, are saved in the cactus folder.
Downstream processing of each single VG archive has been performed on [eddie](https://www.ed.ac.uk/information-services/research-support/research-computing/ecdf/high-performance-computing), the University of Edinburgh high performance computing platform.

## Generation of the indexes on SGE
To generate the indexes on an SGE cluster,you need nodes with at least 500Gb of RAM 
and a minimum of viable space for 2Tb (better 4).
Then, simply run: 

    ./submit.sh

The final graphs will be included into ./GRAPH, including the XG and GCSA indexes.

# Addition of further VCF files
If needed, it is possible to add new variants to a pre-existing graph, [although this is still highly experimental](https://github.com/vgteam/sv-genotyping-paper/issues/6). To do so, proceed as follow:
  1. Create a compliant VCF using the GraphVCF.py script (detailed use of the script can be seen in the submitted script GenerateGraphVCF.sh).
  2. List the new vcf in a file.
  3. Run the AddToGraph.sh script in scripts folder, providing the list of vg graph to expand (-g), the list of VCF to use (-v) and specifying the reference genome to expand (-s)

## Whole genome sequencing analyses
Downstream analysis have been performed using the wrapper [bagpipe](https://bitbucket.org/renzo_tale/bagpipe/src/master/), which allow to perform WGS, ATAC-seq, RRBS and RNA-seq analysis on a SGE/UGE cluster environment. This pipeline also includes script to perform graph genome alignment and variant calling when a graph genome is provided.

