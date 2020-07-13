#!/usr/bin/env nextflow
 
/*
 * Perform liftover from source to target using lastz (different species) 
 * or blat (same species).
 * Specify alignments options in the field alignerparam.
 * Outputs will be saved in the folder specified in outdir.
 * If an annotation is provided as bed/gff, liftover will lift it. 
 * If no annotation is provided, set to 'NO_FILE' or ''
 * Karyotype can be specified as ranges (e.g. 1-22), single different 
 * chromosomes can be added after comma (e.g. 1-22,X,Y,Mt).
 */
params.fai = ''
params.samples = 'samples.txt'
params.genomes = 'genomes.txt'
params.validationFld = '/PATH/TO'
params.ntests = 10000
params.overlap = 0.9


Channel.fromPath(file(params.genomes))
    .splitCsv(header: ['genome'])
    .map{ row-> tuple(row.genome) }
    .into{genomes_ch; genomes_ch2}


Channel.fromPath(file(params.samples))
    .splitCsv(header: ['sampleID', "fullbed", "validbed"])
    .map{ row-> tuple(row.sampleID, row.fullbed, row.validbed) }
    .into{bedfile_ch; bedfile_ch2; bedfile_ch3}


process bootstrapRegion{
    tag "bsRegion"

    input:
    tuple sample, fullbed, validbed from bedfile_ch
    each x from 1..params.ntests

    output:
    tuple sample, x, "myregion.${sample}.${x}.bed" into tests_ch


    script:
    """
    #!/usr/bin/env Rscript
    options(warn=-1, message = FALSE, readr.num_columns = 0, stringsAsFactors = F, scipen=999)
    suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
    suppressPackageStartupMessages(library(GenomicRanges, quietly = TRUE))
    suppressPackageStartupMessages(library(ggpubr, quietly = TRUE))


    myfaifile = read.table("${params.fai}", h=F)
    myfaifile = myfaifile[order(myfaifile[,2]),]
    # Get the number of observations to process
    bedfile = read.table("${fullbed}", h=F)
    positive = nrow(read.table("${validbed}", h=F))
    numPerm = as.numeric("${params.ntest}")
    lengths = as.numeric( abs( bedfile[, 3] - bedfile[, 2] ) )


    getWinds<-function(thisChrom, thisWidth, chr, chr.sizes)
    {
        sample(chr.sizes[chr==thisChrom]-thisWidth,1)
    }

    generateRandomPos <- function(chr, chr.sizes, widths)
    {
        idx = chr.sizes > min(widths)
        random_chr <- sample(x=chr[idx],size=length(widths),prob=chr.sizes[idx],replace=T)
        random_pos <- mapply(getWinds, random_chr, widths, MoreArgs = list(chr[idx], chr.sizes[idx]))
        res <- cbind(random_chr, random_pos, random_pos+widths)
        return(res)
    }

    regions = generateRandomPos(myfaifile[, 'V1'], myfaifile[,'V2'], lengths)
    regions = as.data.frame(cbind(regions, paste("REGION", c(1: nrow(regions)), sep = "")))
    regions[,2] = as.numeric(regions[,2])
    regions[,3] = as.numeric(regions[,3])
    regions[regions[,2] == 0, 2] = 0
    for (i in c(1:nrow(regions))){
        ctgName = regions[i, 1]
        ctgSize = myfaifile[ myfaifile[,1] == ctgName , 2]
        if ( regions[i, 2] > ctgSize ){
            regions[i, 2] = ctgSize
        }
    }
    write.table(regions, "myregion.${sample}.${x}.bed", sep = "\t", row.names = F, col.names = F, quote = F)
    """
}

process intersectRegion{
    tag "intersectRegion"

    input:
    tuple sample, x, regions from tests_ch
    
    output:
    tuple sample, x, "intersected.${sample}.${x}.bed" into intersects_ch


    script:
    """
    bname=\$(basename -s '.bed' ${regions})
    ExtractBed ${regions} \${bname}
    SplitByGenome \${bname}.narrowRegion.bed \$bname.narrow
    while read genome; do
        bedtools intersect -a \${bname}.narrow.\${genome}.bed -b ${params.validationFld}/${sample}/\${genome}.bed -f ${params.overlap} -wa 
    done < ${params.genomes}| sort | uniq > intersected.${sample}.${x}.bed
    """
}


intersects_ch
    .groupTuple(by: [0])
    .set{ grouped_ch }


process getNoverlaps{
    tag "getOverlaps"
    publishDir "${params.outfolder}/LENGTHS/${sample}", mode: 'copy', overwrite: true

    input:
    tuple sample, x, file(intersects) from grouped_ch
    
    output:
    tuple sample, "nintersects.${sample}.txt" into lengths_ch

    script:
    """
    for intersect in ${intersects}; do
        cat \$intersect | sort | uniq | wc -l 
    done > nintersects.${sample}.txt
    """
}

process getPvalues{
    tag "getPval"
    publishDir "${params.outfolder}/PVALUE/${sample}", mode: 'copy', overwrite: true

    input:
    tuple sample, nintersects from lengths_ch
    tuple sample, fullbed, validbed from bedfile_ch2
    
    output:
    tuple sample, "lengths.${sample}.txt" into pvalues_ch

    script:
    """
    #!/usr/bin/env Rscript
    options(warn=-1, message = FALSE, readr.num_columns = 0, stringsAsFactors = F, scipen=999)
    suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
    suppressPackageStartupMessages(library(GenomicRanges, quietly = TRUE))
    suppressPackageStartupMessages(library(ggplot2, quietly = TRUE))
    suppressPackageStartupMessages(library(ggpubr, quietly = TRUE))

    positive = nrow(read.table("${validbed}", h = F))
    found = read.table("${nintersects}", h=F)

    # Calculate permutation values
    Perm_mean = mean(found)
    Perm_sd = sd(found)
    z=(positive-Perm_mean)/Perm_sd 
    z_P<-2*pnorm(-abs(z))

    tot = nrow(regions)
    negative = tot-positive
    pvals = c()
    for (ctrl_pos in found){
        ctrl_neg = tot - ctrl_pos
        csq = chisq.test(c(positive, negative), p=c(ctrl_pos/tot, ctrl_neg/tot))
        pvals = c(pvals, csq[,'p.value'])
    }
    write.table(pvals,"pvals.${sample}.txt", sep = "\t", quote =F, col.names= F, row.names = F)
    final_p = write.table(z_P, "pvalue.${sample}.txt")
    """
}